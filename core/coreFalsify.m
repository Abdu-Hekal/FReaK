function [model,trainset] = coreFalsify(model)

runtime=tic;
%initialization and init assertions
[model, trainset] = initialize(model);
%generate random initial point and inputs
[x0,u] = getRandomXU(model);

falsified = false;
trainIter = 0;
while trainIter < model.maxTrainSize && falsified==false
    [model, trainset] = symbolicRFF(model, trainset, x0, u);
    critX=model.soln.x; critU=model.soln.u;

    %if random trajectory setting selected or critical trajectory is
    %repeated, train with random xu
    if model.trainRand || checkRepeatedTraj(critX,critU)
        [x0,u] = getRandomXU(model);
    else
        x0=critX(1,:)';
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
            robustness = computeRobustness(model.spec(j,1).set,critX,vpa(linspace(0,model.T,size(critX,1)')))
            model.specSolns(model.spec(j,1)).realRob=robustness; %store real robustness value
            check = ~checkStl(model.spec(j,1).set,critX,vpa(linspace(0,model.T,size(critX,1)')));
        end
        if check
            falsified = true;
        end
    end
    trainIter=trainIter+1;
    disp(['iteration completed: ',num2str(trainIter)])
end
%close simulink model
close_system;
%assign solution result
model.soln.falsified=falsified;
model.soln.trainIter=trainIter;
model.soln.runtime=toc(runtime); %record runtime
end

function [x0,u] = getRandomXU(model)
%generate random initial set
x0 = randPoint(model.R0);
%generate random input if model has input.
if ~isempty(model.U)
    all_steps = model.T/model.dt;
    u = randPoint(model.U,all_steps)';
    if model.pulseInput
        u = u.*model.cpBool;
    else %piecewise constant input
        for k=1:length(model.cp)
            uk = model.cp(k);
            u(:,k) = repelem(u(1:uk,k),length(u(:,k))/uk);
        end
    end
    u = [linspace(0,model.T-model.dt,all_steps)',u];

else
    u = [];
end
end

function repeatedTraj = checkRepeatedTraj(critX,critU)
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

%create empty dict to store prev soln (and for each spec)
model.soln=struct;
model.specSolns = dictionary(model.spec,struct);
 %empty cells to store states, inputs and times for training trajectories
trainset.X = {}; trainset.XU={}; trainset.t = {};
end

