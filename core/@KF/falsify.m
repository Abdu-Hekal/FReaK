function [soln,trainset] = falsify(obj)

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
[obj,trainset,soln,specSolns] = initialize(obj);
trainIter = 0;
robustness=inf;

while soln.sims <= obj.maxSims && robustness>=0
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
            specSolns(spec).koopMilp.offsetMap=containers.Map('KeyType', 'double', 'ValueType', 'double');
        end
    end
    % empty trainset at reset and if nonrandom training technique is used, empty trainset after first iter because it is random trajectory.
    if  trainIter ==0 || (trainIter == 1 && obj.trainRand == 0 && obj.rmRand)
        trainset.X = {}; trainset.XU={}; trainset.t = {};
    end

    %if first iter, random trajectory setting selected, or critical trajectory is
    %repeated, retrain with random xu else retrain with prev traj
    if trainIter==0 || obj.trainRand>=2 || ( obj.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(obj,critX,critU, trainset)
        if trainIter>0 && obj.solver.opts.usex0==1 && checkRepeatedTraj(obj,critX,critU, trainset)
            obj.solver.opts.usex0=0; %turn off warmstarting if repeated trajectory returned by solver
            disp('Turned off warmstarting due to repeated solutions')
        end
        [x0,u] = getSampleXU(obj);
        [t, x, simTime] = simulate(obj, x0, u);
        soln.sims = soln.sims+1;
        soln.simTime = soln.simTime+simTime;
        %check if random input falsifies system, and break if it does
        [falsified,robustness]=checkFalsification(x,u,t,obj.spec,obj.inputInterpolation);

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
    %add trajectory to koopman trainset
    trainset.t{end+1} = tak;
    trainset.X{end+1} = xak';
    trainset.XU{end+1} = u(:,2:end)';

    %run autokoopman and learn linearized model
    [koopModel,koopTime] = learnKoopModel(obj, trainset);
    soln.koopTime = soln.koopTime+koopTime;
    % compute reachable set for Koopman linearized model (if reachability is used)
    if obj.reach.on
        reachTime=tic;
        R = reachKoopman(obj,koopModel);
        soln.reachTime=soln.reachTime+toc(reachTime); %time for reachability computation
    else
        R=[];
    end
    % determine most critical reachable set and specification
    try
        optimTime=tic;
        specSolns = critAlpha(obj,R,koopModel,specSolns);
        soln.optimTime=soln.optimTime+toc(optimTime);
    catch
        disp("error encountered whilst setup/solving, resetting training data")
        trainIter=0;
        continue;
    end

    %get critical spec with minimum robustness and corresponding soln struct
    [~,minIndex]=min(specSolns.values.rob);
    keys = specSolns.keys('cell');
    spec=keys{minIndex};
    curSoln=specSolns(spec);

    if curSoln.rob<inf %found a viable solution
        offsetIter = 0;
        while offsetIter <= max(obj.offsetStrat,0) && robustness>=0 %offset once if not falsified
            [critX0, critU] = falsifyingTrajectory(obj,curSoln);
            % run most critical inputs on the real system
            [t, critX, simTime] = simulate(obj, critX0, critU);
            soln.sims = soln.sims+1;
            soln.simTime = soln.simTime+simTime;

            [falsified,robustness,Bdata]=checkFalsification(critX,critU,t,obj.spec,obj.inputInterpolation);

            if abs(robustness)~=inf %there exist a value for robustness for which we can offset
                if obj.trainRand==2 && robustness >= obj.bestSoln.rob %neighborhood training mode
                    %remove last entry because it is not improving the obj
                    trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
                end
                %                 obj=storeBestSoln(obj,robustness,critX,critU,'kf optimization'); %store the best soln so far
                gap = getMilpGap(obj.solver.opts);
                if robustness > gap && abs(obj.offsetStrat) %not falsifed yet, robustness is greater than gap termination criteria for milp solver and an offset mode selected by user. Note that if robustness is less than gap, offset most likely is not benefecial.
                    offsetMap=bReachCulprit(Bdata,spec.set); %get predicates responsible for robustness value
                    if offsetMap.Count > 0 %if there there exists predicates that are culprit for (+ve) robustness
                        obj.solver.opts.usex0=0; %avoid warmstarting if offsetting
                        Sys=specSolns(spec).koopMilp;
                        if offsetIter==0 %if first offset iteration, re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1
                            Sys.offsetMap = offsetMap;
                            if obj.offsetStrat == 1 %if offset strategy in this iteration selected
                                set = spec.set;
                                if ~isequal(set,Sys.stl) || ~obj.useOptimizer
                                    Sys.stl = set;
                                    Sys=setupStl(Sys,~obj.useOptimizer); %encode stl using milp
                                end
                                Sys=optimize(Sys,obj.solver.opts);
                                soln.alpha = value(Sys.alpha); %new alpha value after offset
                            else %obj.offsetStrat == -1: offset next iteration
                                specSolns(spec).koopMilp=Sys;
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
soln.falsified=falsified;
soln.t=t;
soln.x=critX;
soln.u = critU;
soln.runtime=toc(runtime); %record runtime

LogicalStr = {'No', 'Yes'};
vprintf(obj.verb,1,'<-------------------------------------------------------> \n')
vprintf(obj.verb,1,"Falsified: %s \n",LogicalStr{soln.falsified+1})
vprintf(obj.verb,1,"number of simulations %d \n",soln.sims)
vprintf(obj.verb,1,"Time taken %.2f seconds\n",soln.runtime)
end

function [falsified,robustness,Bdata]=checkFalsification(x,u,t,specs,inputInterpolation)
falsified=false;
robustness=inf;
Bdata=NaN;
for ii=1:numel(specs)
    spec=specs(ii);
    % different types of specifications
    if strcmp(spec.type,'unsafeSet')
        falsified = any(spec.set.contains(x'));
    elseif strcmp(spec.type,'safeSet')
        falsified = ~all(spec.set.contains(x')); %check this
    elseif strcmp(spec.type,'logic')
        if ~isempty(u)
            interpU = interp1(u(:,1),u(:,2:end),t,inputInterpolation); %interpolate input at same time points as trajectory
        else
            interpU=u;
        end
        [Bdata,~,robustness] = bReachRob(spec,t,x,interpU');
        if robustness < 0
            falsified=true;
        end
    end
    if falsified
        break;
    end
end
end

function repeatedTraj = checkRepeatedTraj(obj,critX,critU, trainset)
%check if critical initial set & input are the same as found before
repeatedTraj = false;
for r = 1:length(trainset.X)
    if isequal(critX(1,:)',trainset.X{r}(:,1)) && isequal(critU(:,2:end),trainset.XU{r}(:,1:end-1)')
        repeatedTraj = true;
        vprintf(obj.verb,3,"repeated critical trajectory, generating a new random trajectory \n")
        break
    end
end
end

function obj=storeBestSoln(obj,robustness,critX,critU,method)
if robustness < obj.bestSoln.rob
    vprintf(obj.verb,2,"new best robustness!: %.3f after %d simulations due to: %s \n",robustness,soln.sims,method)
    obj.bestSoln.rob = robustness;
    obj.bestSoln.method=method;
    obj.bestSoln.x=critX;
    obj.bestSoln.u=critU;
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

