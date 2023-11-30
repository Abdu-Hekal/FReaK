function [obj,trainset] = falsify(obj)

% falsify - Given a model and a set of specs (safe/unsafe set/stl),
%   perform the core falsification procedure to find a falsifying trajectory
%
% Syntax:
%    [obj, trainset] = falsify(obj)
%
% Description:
%    This function is the core of the falsification procedure, which iteratively
%    searches for a falsifying trajectory by simulating the system, training
%    the Koopman model, and optimizing for critical trajectories. The process
%    continues until a falsifying trajectory is found or a specified timeout is reached.
%
% Inputs:
%    obj - KF object containing the Koopman model and various parameters
%              needed for the falsification process.
%
% Outputs:
%    obj - KF object containing the results of the falsification process,
%              including the falsifying trajectory, simulation time, and other information.
%
%    trainset - Structure containing the training data, including the trajectories and
%               inputs used during the falsification process.
%
% See also:
%
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


runtime=tic;
%initialization and init assertions
[obj, trainset] = initialize(obj);

falsified = false;
trainIter = 0;

while obj.soln.sims <= obj.maxSims && falsified==false
    %timeout
    if toc(runtime) > obj.timeout
        break
    end
    %reset after size of trainset==nResets;
    if numel(trainset.X) == obj.nResets
        trainIter = 0;
        %reset offsets
        for ii=1:numel(obj.spec)
            spec=obj.spec(ii);
            obj.specSolns(spec).lti.offsetMap=containers.Map('KeyType', 'double', 'ValueType', 'double');
        end
    end
    % empty trainset at reset and if nonrandom training technique is used, empty trainset after first iter because it is random trajectory.
    if  trainIter ==0 || (trainIter == 1 && obj.trainRand == 0 && obj.rmRand)
        trainset.X = {}; trainset.XU={}; trainset.t = {};
    end

    %if first iter, random trajectory setting selected, or critical trajectory is
    %repeated, retrain with random xu else retrain with prev traj
    if trainIter==0 || obj.trainRand>=2 || ( obj.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU, trainset)
        if trainIter>0 && obj.solver.opts.usex0==1 && checkRepeatedTraj(critX,critU, trainset)
            obj.solver.opts.usex0=0; %turn off warmstarting if repeated trajectory returned by solver
            disp('Turned off warmstarting due to repeated solutions')
        end
        [x0,u] = getSampleXU(obj);
        tsim = (0:obj.dt:obj.T)'; %define time points for interpolating input
        if ~isempty(u)
            usim = interp1(u(:,1),u(:,2:end),tsim,obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
            usim =  max(obj.U.inf',min(obj.U.sup',usim)); %ensure that extrapolation is within input bounds
            usim = [tsim,usim];
        else
            usim=u; %no input for the model
        end
        [t, x, obj] = simulate(obj, x0, usim);
        %check if random input falsifies system, and break if it does
        falsified = checkFirst(obj,x,usim,t);

        if falsified
            critX=x;
            critU=u;
            break
        end
        tak = (0:obj.ak.dt:obj.T)'; %define autokoopman time points
        xak = interp1(t,x,tak,obj.trajInterpolation); %define autokoopman trajectory points
    else
        xak=interp1(t,critX,tak,obj.trajInterpolation); %pass x0 as full x to avoid simulation again
        u=critU;
    end
    trainset=appendToTrainset(trainset,tak,xak,u);

    %run autokoopman and learn linearized model
    [obj, koopModel] = learnKoopModel(obj, trainset);
    % compute reachable set for Koopman linearized model (if reachability is used)
    if obj.reach.on
        reachTime=tic;
        R = reachKoopman(obj,koopModel);
        obj.soln.reachTime=obj.soln.reachTime+toc(reachTime); %time for reachability computation
    else
        R=[];
    end
    % determine most critical reachable set and specification
    try
        obj = critAlpha(obj,R,koopModel);
    catch
        disp("error encountered whilst setup/solving, resetting training data")
        trainIter=0;
        continue;
    end

    if obj.soln.rob<inf %found a viable solution
        offsetIter = 0;
        while offsetIter <= max(obj.offsetStrat,0) && falsified==false %offset once if not falsified
            [critX0, critU] = falsifyingTrajectory(obj);
            %extrap and interp input for all timesteps T/dt
            if size(critU,2) > 1 %not just time
                usim = interp1(critU(:,1),critU(:,2:end),tsim,obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
                usim =  max(obj.U.inf',min(obj.U.sup',usim)); %ensure that extrapolation is within input bounds
                usim = [tsim,usim];
            else
                usim=critU; %no input for the model
            end
            % run most critical inputs on the real system
            [t, critX, obj] = simulate(obj, critX0, usim);

            %             corePlotFalsify(critU,critX,t,tak,x0,A,B,g,R); %test plot: delete me
            %             corePlotReach(critU,critX,t,tak,x0,A,B,g,R); %test plot: delete me

            if ~isempty(usim)
                interpU = interp1(usim(:,1),usim(:,2:end),t,obj.inputInterpolation); %interpolate input at same time points as trajectory
            else
                interpU=usim;
            end
            spec=obj.soln.spec; %critical spec found with best value of robustness
            % different types of specifications
            if strcmp(spec.type,'unsafeSet')
                falsified = any(spec.set.contains(critX'));
            elseif strcmp(spec.type,'safeSet')
                falsified = ~all(spec.set.contains(critX')); %check this
            elseif strcmp(spec.type,'logic')
                [Bdata,phi,robustness] = bReachRob(spec,t,critX,interpU');

                obj.specSolns(spec).realRob=robustness; %store real robustness value
                falsified = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary
                if obj.trainRand==2 %neighborhood training mode
                    [obj,trainset] = neighborhoodTrain(obj,trainset,robustness,critX,critU);
                end
                gap = getMilpGap(obj.solver.opts);
                if robustness > gap && abs(obj.offsetStrat) %not falsifed yet, robustness is greater than gap termination criteria for milp solver and an offset mode selected by user. Note that if robustness is less than gap, offset most likely is not benefecial.
                    offsetMap=bReachCulprit(Bdata,spec.set); %get predicates responsible for robustness value
                    if offsetMap.Count > 0 %if there there exists predicates that are culprit for (+ve) robustness
                        obj.solver.opts.usex0=0; %avoid warmstarting if offsetting
                        Sys=obj.specSolns(spec).lti;
                        if offsetIter==0 %if first offset iteration, re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1
                            Sys.offsetMap = offsetMap;
                            if obj.offsetStrat == 1 %if offset strategy in this iteration selected
                                set = spec.set;
                                if ~isequal(set,Sys.stl) || ~obj.useOptimizer
                                    Sys.stl = set;
                                    Sys=setupStl(Sys,~obj.useOptimizer); %encode stl using milp
                                end
                                Sys=optimize(Sys,obj.solver.opts);
                                obj.soln.alpha = value(Sys.alpha); %new alpha value after offset
                            else %obj.offsetStrat == -1: offset next iteration
                                obj.specSolns(spec).lti=Sys;
                            end
                            % TODO: if offset gives better val of robustness, should we pass
                            %                     % it as training data instead? should we pass both?
                        end
                    else
                        break
                    end
                else
                    break
                end
            end
            offsetIter = offsetIter+1;
        end
    end
    %clear previous solution from yalmip (makes warmstart feasible)
    if obj.solver.opts.usex0
        yalmip('clearsolution')
    end
    trainIter=trainIter+1;
end
%close simulink obj
close_system('ErrorIfShadowed',0);
%assign solution result
obj.soln.t=t;
obj.soln.x=critX;
obj.soln.u = critU;
obj.soln.falsified=falsified;
obj.soln.runtime=toc(runtime); %record runtime
end

function repeatedTraj = checkRepeatedTraj(critX,critU, trainset)
%check if critical initial set & input are the same as found before
repeatedTraj = false;
for r = 1:length(trainset.X)
    if isequal(critX(1,:)',trainset.X{r}(:,1)) && isequal(critU(:,2:end),trainset.XU{r}(:,1:end-1)')
        repeatedTraj = true;
        disp("repeated critical trajectory, generating a new random trajectory")
        break
    end
end
end

function falsified = checkFirst(obj,x,usim,t)
if ~isempty(usim)
    interpU = interp1(usim(:,1),usim(:,2:end),t,obj.inputInterpolation); %interpolate input at same time points as trajectory
else
    interpU=usim;
end
for ii=1:numel(obj.spec)
    spec=obj.spec(ii);
    % different types of specifications
    if strcmp(spec.type,'unsafeSet')
        falsified = any(spec.set.contains(x'));
    elseif strcmp(spec.type,'safeSet')
        falsified = ~all(spec.set.contains(x')); %check this
    elseif strcmp(spec.type,'logic')
        [~,~,robustness] = bReachRob(spec,t,x,interpU');
        falsified = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary
    end
end
end

function trainset=appendToTrainset(trainset,t,x,u)
%add trajectory to koopman trainset
trainset.t{end+1} = t;
trainset.X{end+1} = x';
trainset.XU{end+1} = u(:,2:end)';
end

function [obj,trainset] = neighborhoodTrain(obj,trainset,robustness,critX,critU)
%training with best neighborhood trajectories
if robustness < obj.bestSoln.rob
    disp("new best robustness found!")
    disp(robustness)
    obj.bestSoln.rob = robustness;
    obj.bestSoln.x=critX;
    obj.bestSoln.u=critU;
else
    %remove last entry because it is not improving the obj
    trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
end
end

function gap=getMilpGap(opts)
if opts.solver == "gurobi"
    gap=opts.gurobi.MIPGapAbs;
elseif opts.solver == "cplex"
    gap=opts.cplex.epagap;
else
    gap=0.1;
end
end

