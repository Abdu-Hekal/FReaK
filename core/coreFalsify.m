function [falsified, trainset, crit_x, train_iter] = coreFalsify(model, max_train_size)

%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');
%initialization and init assertions
model = initialize(model);
%generate random initial point and inputs
[x0,u] = gen_random_xu(model);

trainset.X = {}; trainset.XU={}; trainset.t = {}; %empty cells to store states, inputs and times for training trajectories

falsified = false;
train_iter = 0;
while train_iter < max_train_size && falsified==false
    [trainset, crit_x, crit_u] = symbolicRFF(model, trainset, x0, u);
    %retrain with initial set & input as the critical values found in prev iteration.
    %If the critical values are the same as previous, generate new random
    %values
    repeated_traj = false;
    for r = 1:length(trainset.X)
        if isequal(crit_x(1,:)',trainset.X{r}(:,1)) && isequal(crit_u(:,2:end),trainset.XU{r}(:,1:end-1)')
            repeated_traj = true;
            break
        end
    end
    if true %repeated_traj
        disp("repeated critical trajectory, generating a new random trajectory")
        [x0,u] = gen_random_xu(model);
    else
        x0=crit_x(1,:)';
        u=crit_u;
    end
    %     disp(crit_u)

    for j = 1:size(model.spec,1)
        % different types of specifications
        if strcmp(model.spec(j,1).type,'unsafeSet')
            check = any(model.spec(j,1).set.contains(crit_x'));
        elseif strcmp(model.spec(j,1).type,'logic')
            robustness = computeRobustness(model.spec(j,1).set,crit_x,vpa(linspace(0,model.T,size(crit_x,1)')))
            check = ~checkStl(model.spec(j,1).set,crit_x,vpa(linspace(0,model.T,size(crit_x,1)')));
        end
        if check
            falsified = true;
        end
    end
    train_iter=train_iter+1;
    disp(['iteration completed: ',num2str(train_iter)])
end
%close simulink model
close_system;
end

function [x0,u] = gen_random_xu(model)
%generate random initial set
x0 = randPoint(model.R0);
%generate random input if model has input.
if ~isempty(model.U)
    all_steps = model.T/model.dt;
    u = randPoint(model.U,all_steps)';
    if model.pulse_input
        u = u.*model.cp_bool;
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

function model = initialize(model)
assert(isa(model.sim, 'string') | isa(model.sim,"char")| isa(model.sim,'function_handle'), 'model.sim must be a (1)string: name of simulink model or a (2)function handle')
assert(isa(model.R0, 'interval'), 'Initial set (model.R0) must be defined as an CORA interval')
assert(~isempty(model.T) & isnumeric(model.T), 'Time horizon (model.T) must be defined as a numeric')
assert(~isempty(model.dt) & isnumeric(model.dt), 'Time step (model.dt) must be defined as a numeric')
assert(isa(model.spec, 'specification'), 'Falsifying spec (model.spec) must be defined as a CORA specification')

if ~isempty(model.U) %check if simulink model has inputs
    assert(isa(model.U, 'interval'), 'Input (model.U) must be defined as an CORA interval')
    %if no control points defined, set as control point at every step dt
    if isempty(model.cp)
        model.cp=model.T/model.dt*ones(1,length(model.U));
    end

    assert(length(model.U)==length(model.cp),'Number of control points (model.cp) must be equal to number of inputs (model.U)')

    all_steps = model.T/model.dt;
    assert(floor(all_steps)==all_steps,'Time step (dt) must be a factor of Time horizon (T)')

    model.cp_bool = zeros(all_steps,length(model.U));
    for k=1:length(model.cp)
        step = (model.T/model.dt)/model.cp(k);
        assert(floor(step)==step,'number of control points (cp) must be a factor of T/dt')
        model.cp_bool(1:step:end,k) = 1;
    end
end

%create empty dict to store soln for each spec
model.spec_soln = dictionary(model.spec,{{}});
end

