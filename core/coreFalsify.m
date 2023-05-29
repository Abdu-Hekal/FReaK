function [kfModel,trainset] = coreFalsify(kfModel)

runtime=tic;
%initialization and init assertions
[kfModel, trainset] = initialize(kfModel);
%generate random initial point and inputs
[x0,u] = getSampleXU(kfModel);

%test: pretrain with n traj
% for i=1:400
%     [x0,u] = getSampleXU(kfModel);
%     [trainset.t{end+1}, x, kfModel] = simulate(kfModel, x0, u);
%     trainset.X{end+1} = x';
%     trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];
% end

falsified = false;
trainIter = 0;
while trainIter <= kfModel.maxTrainSize && falsified==false

    %run autokoopman and learn linearized model
    [kfModel, trainset, A, B, g] = symbolicRFF(kfModel, trainset, x0, u);
    % finf predicted falsifying initial set and inputs
    % compute reachable set for Koopman linearized model
    R = reachKoopman(A,B,g,kfModel);
    % determine most critical reachable set and specification
    kfModel = critAlpha(R,kfModel);

    refineIter = 0;
    while refineIter <= 1 && falsified==false %refine once if not falsified
        [critX0, critU] = falsifyingTrajectory(kfModel);
        % run most critical inputs on the real system
        [t, critX, kfModel] = simulate(kfModel, critX0, critU);

        spec=kfModel.soln.spec; %critical spec found with best value of robustness
        % different types of specifications
        if strcmp(spec.type,'unsafeSet')
            falsified = any(spec.set.contains(critX'));
        elseif strcmp(spec.type,'safeSet')
            falsified = ~all(spec.set.contains(critX')); %check this
        elseif strcmp(spec.type,'logic')
            %test plot: delete me
            plot(critX(1:400,1),critX(1:400,2),'g','LineWidth',2)
            drawnow;
            robustness = computeRobustness(spec.set,critX,vpa(linspace(0,kfModel.T,size(critX,1)')))
            %  breachRob = bReachRob(kfModel.spec,kfModel.soln.x,kfModel.soln.t)
            kfModel.specSolns(spec).realRob=robustness; %store real robustness value
            falsified = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary

            if ~falsified && refineIter < 1 %not falsifed, robustness > 0: re-solve with offset
                Sys=kfModel.specSolns(spec).lti;
                Sys.offset = robustness;
                Sys.offsetCount = getRobOffset(spec.set,critX,vpa(linspace(0,kfModel.T,size(critX,1)')),robustness);
                Sys=optimize(Sys);
                kfModel.soln.alpha = value(Sys.Alpha); %new alpha value after offset
            end
            refineIter = refineIter+1;

            %training with best neighborhood trajectories
            if kfModel.trainRand==2
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
        end
    end
    trainIter=trainIter+1;

    %retrain: if random trajectory setting selected or critical trajectory is
    %repeated, retrain with random xu else retrain with prev traj
    if kfModel.trainRand>=2 || ( kfModel.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU, trainset)
        [x0,u] = getSampleXU(kfModel);
    else
        x0=critX; %pass x0 as full x to avoid simulation again
        u=critU;
    end
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

%add dt to autokoopman settings as well
kfModel.ak.dt=kfModel.dt;
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

    all_steps = kfModel.T/kfModel.dt;
    assert(floor(all_steps)==all_steps,'Time step (dt) must be a factor of Time horizon (T)')

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

