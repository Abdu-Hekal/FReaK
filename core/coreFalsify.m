function [kfModel,trainset] = coreFalsify(kfModel)

runtime=tic;
%initialization and init assertions
[kfModel, trainset] = initialize(kfModel);

falsified = false;
trainIter = 0;
epsilon = eps; %epsilon value to add to offset
while trainIter <= kfModel.maxTrainSize && falsified==false

    %empty trainset after first iter because it is random trajectory
    if trainIter == 1
        trainset.X = {}; trainset.XU={}; trainset.t = {};
    end

    %if first iter, random trajectory setting selected, or critical trajectory is
    %repeated, retrain with random xu else retrain with prev traj
    if trainIter==0 || kfModel.trainRand>=2 || ( kfModel.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU, trainset)
        [x0,u] = getSampleXU(kfModel);
        [t, x, kfModel] = simulate(kfModel, x0, u);
    else
        x=critX; %pass x0 as full x to avoid simulation again
        u=critU;
    end
    abstr = kfModel.ak.dt/kfModel.dt; %define abstraction ratio
    t=t(1:abstr:end); x=x(1:abstr:end,:); u=u(1:abstr:end,:);
    trainset=appendToTrainset(trainset,t,x,u);

    %run autokoopman and learn linearized model
    [kfModel, A, B, g] = symbolicRFF(kfModel, trainset);
    % find predicted falsifying initial set and inputs
    % compute reachable set for Koopman linearized model
    R = reachKoopman(A,B,g,kfModel);
    % determine most critical reachable set and specification
    kfModel = critAlpha(R,kfModel);

    offsetIter = 0;
    while offsetIter <= max(kfModel.offsetStrat,0) && falsified==false %offset once if not falsified
        [critX0, critU] = falsifyingTrajectory(kfModel);
        %repeat input for all timesteps T/dt
        oldU=critU; %test deleteme
        critU=repelem(critU,abstr,1);
        % run most critical inputs on the real system
        [t, critX, kfModel] = simulate(kfModel, critX0, critU);
        testDraw(oldU,critX,t,t(1:abstr:end),x0,A,B,g,R); %test plot: delete me

        spec=kfModel.soln.spec; %critical spec found with best value of robustness
        % different types of specifications
        if strcmp(spec.type,'unsafeSet')
            falsified = any(spec.set.contains(critX'));
        elseif strcmp(spec.type,'safeSet')
            falsified = ~all(spec.set.contains(critX')); %check this
        elseif strcmp(spec.type,'logic')
            [Bdata,phi,robustness] = bReachRob(spec,critX,t);
            kfModel.specSolns(spec).realRob=robustness; %store real robustness value
            falsified = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary
            if kfModel.trainRand==2 && numel(trainset.X)>1 %neighborhood training mode and we have more than 1 trajectory in our trainset
                [kfModel,trainset] = neighborhoodTrain(kfModel,trainset,robustness,critX,critU);
            end

            if ~falsified && abs(kfModel.offsetStrat) %not falsifed yet and an offset mode selected by user
                newOffsetCount = bReachCulprit(Bdata,phi,robustness); %get no. of predicate responsible for robustness value
                if newOffsetCount > 0 %there exists an individual predicate that is culprit for (+ve) robustness
                    if offsetIter==0 %if first offset iteration, re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1
                        Sys=kfModel.specSolns(spec).lti;
                        if Sys.offsetCount == newOffsetCount %if same offset count as b4 offset (i.e. same inequality)
                            Sys.offset =  Sys.offset+robustness+epsilon;
                        else
                            Sys.offset = robustness + epsilon;
                        end
                        Sys.offsetCount = newOffsetCount;
                        if kfModel.offsetStrat == 1 %if offset strategy in this iteration selected
                            if ~kfModel.useOptimizer %if no optimizer object, setup stl with hardcoded offset
                                Sys=setupStl(Sys,true);
                            end
                            Sys=optimize(Sys,kfModel.solver.opts);
                            kfModel.soln.alpha = value(Sys.Alpha); %new alpha value after offset
                        else %kfModel.offsetStrat == -1: offset next iteration
                            kfModel.specSolns(spec).lti=Sys;
                        end
                        % TODO: if offset gives better val of robustness, should we pass
                        %                     % it as training data instead? should we pass both?
                    else %if already offset and failed to falsify, increase epsilon
                        %check first if same offset count as b4 offset (i.e. same inequality)
                        if Sys.offsetCount == newOffsetCount
                            epsilon = epsilon + robustness;
                        else
                            epsilon = eps; %reset epsilon cause we now violate a different inequality
                        end
                    end
                end
                offsetIter = offsetIter+1;
            end
        end
    end
    trainIter=trainIter+1;
end
%close simulink kfModel
close_system;
%assign solution result
kfModel.soln.t=t;
kfModel.soln.x=critX;
kfModel.soln.u = critU;
kfModel.soln.falsified=falsified;
kfModel.soln.trainIter=trainIter;
kfModel.soln.sims=kfModel.soln.sims-1; %last sim that finds critical trace not counted?
kfModel.soln.runtime=toc(runtime); %record runtime
fprintf("number of simulations to falsify %d \n",kfModel.soln.sims)
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
if isempty(kfModel.ak.dt)
    kfModel.ak.dt=kfModel.dt;
else
    abstr = kfModel.ak.dt/kfModel.dt; %define abstraction ratio
    assert(floor(abstr)==abstr,'Time step of koopman (ak.dt) must be a multiple of dt')
    allAbstrSteps = kfModel.T/kfModel.ak.dt;
    assert(floor(allAbstrSteps)==allAbstrSteps,'Time step of koopman (ak.dt) must be a factor of Time horizon (T)')
end

%set solver timestep if it is not set, else check it is compliant.
if isempty(kfModel.solver.dt)
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

% initialize seed
rng(0)

if ~isempty(kfModel.U) %check if kfModel has inputs
    assert(isa(kfModel.U, 'interval'), 'Input (kfModel.U) must be defined as an CORA interval')
    %if no control points defined, set as control point at every step dt
    if isempty(kfModel.cp)
        kfModel.cp=kfModel.T/kfModel.dt*ones(1,length(kfModel.U));
    end

    assert(length(kfModel.U)==length(kfModel.cp),'Number of control points (kfModel.cp) must be equal to number of inputs (kfModel.U)')

    kfModel.cpBool = zeros(all_steps,length(kfModel.U));
    for k=1:length(kfModel.cp)
        step = (kfModel.T/kfModel.dt)/kfModel.cp(k);
        assert(floor(step)==step,'number of control points (cp) must be a factor of T/dt')
        kfModel.cpBool(1:step:end,k) = 1;
    end
end

%create empty struct to store prev soln
kfModel.soln=struct;
kfModel.soln.koopTime=0; kfModel.soln.milpSetupTime=0; kfModel.soln.milpSolvTime=0;
kfModel.soln.sims=0;
%create empty struct to store best soln
kfModel.bestSoln=struct; kfModel.bestSoln.rob=inf;
%create empty dict to store prev soln for each spec
kfModel.specSolns = dictionary(kfModel.spec,struct);
%empty cells to store states, inputs and times for training trajectories
trainset.X = {}; trainset.XU={}; trainset.t = {};
end

function trainset=appendToTrainset(trainset,t,x,u)
%add trajectory to koopman trainset
trainset.t{end+1} = t;
trainset.X{end+1} = x';
%append zeros at end to account for last time point (which has no
%inputs), but length must be consistant with trajectory states
trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];
end

function [kfModel,trainset] = neighborhoodTrain(kfModel,trainset,robustness,critX,critU)
%training with best neighborhood trajectories
if robustness < kfModel.bestSoln.rob
    kfModel.bestSoln.rob = robustness;
    kfModel.bestSoln.x=critX;
    kfModel.bestSoln.u=critU;
else
    %remove last entry because it is not improving the kfModel
    trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
end
disp(kfModel.bestSoln.rob)
end

function testDraw(critU,critX,t,xt,x0,A,B,g,R)
plotVars=[3]; %[3];
drawu=critU(:,2:end)';
x = g(x0);
for i = 1:size(drawu,2)
    x = [x, A*x(:,end) + B*drawu(:,i)];
end
figure; hold on; box on;
if ~any(size(plotVars)>[1,1]) %singular plot var, plot against time
    plot(xt,x(plotVars(1),1:end),'r','LineWidth',2);
    plot(t,critX(1:end,plotVars(1)),'g','LineWidth',2)
else
    for i=1:size(drawu,2)
        plot(R.zono{i},plotVars)
    end
    plot(x(plotVars(1),1:end),x(plotVars(2),1:end),'r','LineWidth',2);
    plot(critX(1:end,plotVars(1)),critX(1:end,plotVars(2)),'g','LineWidth',2)
end
drawnow;
end

