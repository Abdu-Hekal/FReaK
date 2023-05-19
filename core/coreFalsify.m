function [model,trainset] = coreFalsify(model)

runtime=tic;
%initialization and init assertions
[model, trainset] = initialize(model);
%generate random initial point and inputs
[x0,u] = getSampleXU(model);

%test: pretrain with n traj
% for i=1:400
%     [x0,u] = getSampleXU(model);
%     [trainset.t{end+1}, x, model] = simulate(model, x0, u);
%     trainset.X{end+1} = x';
%     trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];
% end

falsified = false;
trainIter = 0;
while trainIter < model.maxTrainSize && falsified==false
    [model, trainset] = symbolicRFF(model, trainset, x0, u);
    critX=model.soln.x; critU=model.soln.u;

    %if random trajectory setting selected or critical trajectory is
    %repeated, train with random xu
    if model.trainRand==2 || ( model.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU, trainset)
        [x0,u] = getSampleXU(model);
    else
        x0=critX; %pass x0 as full x to avoid simulation again
        u=critU;
    end
    %     disp(critU)

    %check for all spec, if by luck trace falsifies a different spec than crit spec ;)
    for j = 1:size(model.spec,1)
        % different types of specifications
        if strcmp(model.spec(j,1).type,'unsafeSet')
            check = any(model.spec(j,1).set.contains(critX'));
        elseif strcmp(model.spec(j,1).type,'safeSet')
            check = ~all(model.spec(j,1).set.contains(critX')); %check this
        elseif strcmp(model.spec(j,1).type,'logic')
            plot(critX(1:400,1),critX(1:400,2),'g','LineWidth',2)
            drawnow;
            robustness = computeRobustness(model.spec(j,1).set,critX,vpa(linspace(0,model.T,size(critX,1)')))
            if robustness < model.bestSoln.rob
                model.bestSoln.rob = robustness;
                model.bestSoln.x=critX;
                model.bestSoln.u=critU;
            else
                %remove last entry because it is not improving the model
                trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
            end
            disp(model.bestSoln.rob)
%             breachRob = bReachRob(model.spec,model.soln.x,model.soln.t)
            model.specSolns(model.spec(j,1)).realRob=robustness; %store real robustness value
            check = ~isreal(sqrt(robustness)); %sqrt of -ve values are imaginary
            %             check= ~isreal(sqrt(breachRob));
        end
        if check
            falsified = true;
        end
    end
    trainIter=trainIter+1;
end
%close simulink model
close_system;
%assign solution result
model.soln.falsified=falsified;
model.soln.trainIter=trainIter;
model.soln.sims=model.soln.sims-1; %last sim that finds critical trace not counted?
model.soln.runtime=toc(runtime); %record runtime
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

function [model, trainset] = initialize(model)
%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

assert(isa(model.sim, 'string') | isa(model.sim,"char")| isa(model.sim,'function_handle'), 'model.sim must be a (1)string: name of simulink model or a (2)function handle')
assert(isa(model.R0, 'interval'), 'Initial set (model.R0) must be defined as an CORA interval')
assert(~isempty(model.T) & isnumeric(model.T), 'Time horizon (model.T) must be defined as a numeric')
assert(~isempty(model.dt) & isnumeric(model.dt), 'Time step (model.dt) must be defined as a numeric')
assert(isa(model.spec, 'specification'), 'Falsifying spec (model.spec) must be defined as a CORA specification')

% clear yalmip
yalmip('clear')

% initialize seed
rng(0)

if ~isempty(model.U) %check if model has inputs
    assert(isa(model.U, 'interval'), 'Input (model.U) must be defined as an CORA interval')
    %if no control points defined, set as control point at every step dt
    if isempty(model.cp)
        model.cp=model.T/model.dt*ones(1,length(model.U));
    end

    assert(length(model.U)==length(model.cp),'Number of control points (model.cp) must be equal to number of inputs (model.U)')

    all_steps = model.T/model.dt;
    assert(floor(all_steps)==all_steps,'Time step (dt) must be a factor of Time horizon (T)')

    model.cpBool = zeros(all_steps,length(model.U));
    for k=1:length(model.cp)
        step = (model.T/model.dt)/model.cp(k);
        assert(floor(step)==step,'number of control points (cp) must be a factor of T/dt')
        model.cpBool(1:step:end,k) = 1;
    end
end

%create empty struct to store prev soln
model.soln=struct;
model.soln.koopTime=0; model.soln.milpSetupTime=0; model.soln.milpSolvTime=0;
model.soln.sims=0;
%create empty struct to store best soln
model.bestSoln=struct; model.bestSoln.rob=inf;
%create empty dict to store prev soln for each spec
model.specSolns = dictionary(model.spec,struct);
%empty cells to store states, inputs and times for training trajectories
trainset.X = {}; trainset.XU={}; trainset.t = {};
end

