function [kfModel,trainset] = coreFalsify(kfModel)

runtime=tic;
%initialization and init assertions
[kfModel, trainset] = initialize(kfModel);

falsified = false;
trainIter = 0;

while kfModel.soln.sims <= kfModel.maxSims && falsified==false
    %timeout
    if toc(runtime) > kfModel.timeout
        break
    end
    %reset after size of trainset==nResets;
    if numel(trainset.X) == kfModel.nResets
        trainIter = 0;
        %reset offsets
        for ii=1:numel(kfModel.spec)
            spec=kfModel.spec(ii);
            kfModel.specSolns(spec).lti.offsetMap=containers.Map('KeyType', 'double', 'ValueType', 'double');
        end
    end
    % empty trainset at reset and if nonrandom training technique is used, empty trainset after first iter because it is random trajectory.
    if  trainIter ==0 || (trainIter == 1 && kfModel.trainRand == 0 && kfModel.rmRand)
        trainset.X = {}; trainset.XU={}; trainset.t = {};
    end

    kfModel=setCpBool(kfModel);
    %if first iter, random trajectory setting selected, or critical trajectory is
    %repeated, retrain with random xu else retrain with prev traj
    if trainIter==0 || kfModel.trainRand>=2 || ( kfModel.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU, trainset)
        if trainIter>0 && kfModel.solver.opts.usex0==1 && checkRepeatedTraj(critX,critU, trainset)
            kfModel.solver.opts.usex0=0; %turn off warmstarting if repeated trajectory returned by solver
            disp('Turned off warmstarting due to repeated solutions')
        end
        [x0,u] = getSampleXU(kfModel);
        tsim = (0:kfModel.dt:kfModel.T)'; %define time points for interpolating input
        if ~isempty(u)
            usim = interp1(u(:,1),u(:,2:end),tsim,kfModel.inputInterpolation,"extrap"); %interpolate and extrapolate input points
            usim =  max(kfModel.U.inf',min(kfModel.U.sup',usim)); %ensure that extrapolation is within input bounds
            usim = [tsim,usim];
        else
            usim=u; %no input for the model
        end
        [t, x, kfModel] = simulate(kfModel, x0, usim);
        %check if random input falsifies system, and break if it does
        falsified = checkFirst(kfModel,x,usim,t);

        if falsified
            critX=x;
            critU=u;
            break
        end
        tak = (0:kfModel.ak.dt:kfModel.T)'; %define autokoopman time points
        xak = interp1(t,x,tak,kfModel.trajInterpolation); %define autokoopman trajectory points
    else
        xak=interp1(t,critX,tak,kfModel.trajInterpolation); %pass x0 as full x to avoid simulation again
        u=critU;
    end
    trainset=appendToTrainset(trainset,tak,xak,u);

    %     try
    %run autokoopman and learn linearized model
    [kfModel, A, B, g] = symbolicRFF(kfModel, trainset);
    % compute reachable set for Koopman linearized model (if reachability is used)
    if kfModel.reach
        R = reachKoopman(A,B,g,kfModel);
    else
        R=[];
    end
    % determine most critical reachable set and specification
    kfModel = critAlpha(R,A,B,g,kfModel);
    %     catch
    %         disp("error encountered whilst setup/solving, resetting training data")
    %         trainIter=0;
    %         continue;
    %     end

    if kfModel.soln.rob<inf %found a viable solution
        offsetIter = 0;
        while offsetIter <= max(kfModel.offsetStrat,0) && falsified==false %offset once if not falsified
            [critX0, critU] = falsifyingTrajectory(kfModel);
            %extrap and interp input for all timesteps T/dt
            if size(critU,2) > 1 %not just time
                usim = interp1(critU(:,1),critU(:,2:end),tsim,kfModel.inputInterpolation,"extrap"); %interpolate and extrapolate input points
                usim =  max(kfModel.U.inf',min(kfModel.U.sup',usim)); %ensure that extrapolation is within input bounds
                usim = [tsim,usim];
            else
                usim=critU; %no input for the model
            end
            % run most critical inputs on the real system
            [t, critX, kfModel] = simulate(kfModel, critX0, usim);

            %             corePlotFalsify(critU,critX,t,tak,x0,A,B,g,R); %test plot: delete me
            %             corePlotReach(critU,critX,t,tak,x0,A,B,g,R); %test plot: delete me

            if ~isempty(usim)
                interpU = interp1(usim(:,1),usim(:,2:end),t,kfModel.inputInterpolation); %interpolate input at same time points as trajectory
            else
                interpU=usim;
            end
            spec=kfModel.soln.spec; %critical spec found with best value of robustness
            % different types of specifications
            if strcmp(spec.type,'unsafeSet')
                falsified = any(spec.set.contains(critX'));
            elseif strcmp(spec.type,'safeSet')
                falsified = ~all(spec.set.contains(critX')); %check this
            elseif strcmp(spec.type,'logic')
                [Bdata,phi,robustness] = bReachRob(spec,t,critX,interpU');

                kfModel.specSolns(spec).realRob=robustness; %store real robustness value
                falsified = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary
                if kfModel.trainRand==2 %neighborhood training mode
                    [kfModel,trainset] = neighborhoodTrain(kfModel,trainset,robustness,critX,critU);
                end
                gap = getMilpGap(kfModel.solver.opts);
                if robustness > gap && abs(kfModel.offsetStrat) %not falsifed yet, robustness is greater than gap termination criteria for milp solver and an offset mode selected by user. Note that if robustness is less than gap, offset most likely is not benefecial.
                    offsetMap=bReachCulprit(Bdata,spec.set); %get predicates responsible for robustness value
                    if offsetMap.Count > 0 %if there there exists predicates that are culprit for (+ve) robustness
                        kfModel.solver.opts.usex0=0; %avoid warmstarting if offsetting
                        Sys=kfModel.specSolns(spec).lti;
                        if offsetIter==0 %if first offset iteration, re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1
                            Sys.offsetMap = offsetMap;
                            if kfModel.offsetStrat == 1 %if offset strategy in this iteration selected
                                set = spec.set;
                                bluStl = coraBlustlConvert(set); %convert from cora syntax to blustl
                                if ~isequal(bluStl,Sys.stl) || ~kfModel.useOptimizer
                                    Sys.stl = bluStl;
                                    Sys=setupStl(Sys,~kfModel.useOptimizer); %encode stl using milp
                                end
                                Sys=optimize(Sys,kfModel.solver.opts);
                                kfModel.soln.alpha = value(Sys.alpha); %new alpha value after offset
                            else %kfModel.offsetStrat == -1: offset next iteration
                                kfModel.specSolns(spec).lti=Sys;
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
    if kfModel.solver.opts.usex0
        yalmip('clearsolution')
    end
    trainIter=trainIter+1;
end
%close simulink kfModel
close_system(' ErrorIfShadowed',0);
%assign solution result
kfModel.soln.t=t;
kfModel.soln.x=critX;
kfModel.soln.u = critU;
kfModel.soln.falsified=falsified;
kfModel.soln.runtime=toc(runtime); %record runtime
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

function [kfModel, trainset] = initialize(kfModel)
%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

assert(isa(kfModel.model, 'string') | isa(kfModel.model,"char")| isa(kfModel.model,'function_handle'), 'kfModel.model must be a (1)string: name of simulink kfModel or a (2)function handle')
assert(isa(kfModel.R0, 'interval'), 'Initial set (kfModel.R0) must be defined as an CORA interval')
assert(~isempty(kfModel.T) & isnumeric(kfModel.T), 'Time horizon (kfModel.T) must be defined as a numeric')
assert(~isempty(kfModel.dt) & isnumeric(kfModel.dt), 'Time step (kfModel.dt) must be defined as a numeric')
assert(isa(kfModel.spec, 'specification'), 'Falsifying spec (kfModel.spec) must be defined as a CORA specification')
all_steps = kfModel.T/kfModel.dt;
assert(floor(all_steps)==all_steps,'Time step (dt) must be a factor of Time horizon (T)')

%set autokoopman timestep if it is not set, else check it is compliant.
if ~isfield(kfModel.ak,'dt')
    kfModel.ak.dt=kfModel.dt;
else
    allAbstrSteps = kfModel.T/kfModel.ak.dt;
    assert(floor(allAbstrSteps)==allAbstrSteps,'Time step of koopman (ak.dt) must be a factor of Time horizon (T)')
end

%set solver timestep if it is not set, else check it is compliant.
if ~isfield(kfModel.solver,'dt')
    kfModel.solver.dt=kfModel.ak.dt;
else
    abstr = kfModel.solver.dt/kfModel.ak.dt; %define abstraction ratio
    assert(floor(abstr)==abstr,'Time step of solver (solver.dt) must be a multiple of koopman time step (ak.dt)')
    allAbstrSteps = kfModel.T/kfModel.solver.dt;
    assert(floor(allAbstrSteps)==allAbstrSteps,'Time step of solver (solver.dt) must be a factor of Time horizon (T)')
end

% ensure that autokoopman rank is an integer
kfModel.ak.rank=int64(kfModel.ak.rank);

% clear yalmip
yalmip('clear')

if ~isempty(kfModel.U) %check if kfModel has inputs
    assert(isa(kfModel.U, 'interval'), 'Input (kfModel.U) must be defined as an CORA interval')
    %if no control points defined, set as control point at every step dt
    if isempty(kfModel.cp)
        kfModel.cp=kfModel.T/kfModel.dt*ones(1,length(kfModel.U));
    end

    assert(length(kfModel.U)==length(kfModel.cp),'Number of control points (kfModel.cp) must be equal to number of inputs (kfModel.U)')
end

%reset struct to store prev soln
kfModel.soln=struct;
kfModel.soln.koopTime=0; kfModel.soln.milpSetupTime=0; kfModel.soln.milpSolvTime=0; kfModel.soln.simTime=0;
kfModel.soln.sims=0;
%reset struct to store best soln
kfModel.bestSoln=struct; kfModel.bestSoln.rob=inf; kfModel.bestSoln.timeRob=inf;
%reset dict to store prev soln for each spec
kfModel.specSolns = dictionary(kfModel.spec,struct);
%empty struct to store training data
trainset.X = {}; trainset.XU={}; trainset.t = {};
end

function falsified = checkFirst(kfModel,x,usim,t)
if ~isempty(usim)
    interpU = interp1(usim(:,1),usim(:,2:end),t,kfModel.inputInterpolation); %interpolate input at same time points as trajectory
else
    interpU=usim;
end
for ii=1:numel(kfModel.spec)
    spec=kfModel.spec(ii);
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

function [kfModel,trainset] = neighborhoodTrain(kfModel,trainset,robustness,critX,critU)
%training with best neighborhood trajectories
if robustness < kfModel.bestSoln.rob
    disp("new best robustness found!")
    disp(robustness)
    kfModel.bestSoln.rob = robustness;
    kfModel.bestSoln.x=critX;
    kfModel.bestSoln.u=critU;
else
    %remove last entry because it is not improving the kfModel
    trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
end
end

function kfModel=setCpBool(kfModel)
all_steps = kfModel.T/kfModel.ak.dt;
kfModel.cpBool = zeros(all_steps,length(kfModel.U));
for k=1:length(kfModel.cp)
    step = (kfModel.T/kfModel.ak.dt)/kfModel.cp(k);
    step = max(1,step); %if step<1, then we have more control points than ak steps, set control points equal to number of steps
    assert(floor(step)==step,'number of control points (cp) must be a factor of T/ak.dt')
    kfModel.cpBool(1:step:end,k) = 1;
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

