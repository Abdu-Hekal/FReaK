function [falsified, trainset, crit_x, train_iter] = coreFalsify(model, max_train_size)

%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

%generate random initial set and input
x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
u = [linspace(0,model.T-model.T/model.N,model.N)',((model.U.sup-model.U.inf).*rand(1,model.N)+model.U.inf)']; %check this

trainset.X = {}; trainset.XU={}; trainset.t = {}; %empty cells to store states, inputs and times for training trajectories

falsified = false;
i = 0;

%headstart with n trajs
% for k=1:20
%     [trainset.t{end+1}, x] = run_simulation(model.name, model.T, x0, u);
%     trainset.X{end+1} = x';
%     %append zeros at end to account for last time point (which has no
%     %inputs), but length must be consistant with trajectory states
%     trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];
%     x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
%     u = [linspace(0,model.T-model.T/model.N,model.N)',((model.U.sup-model.U.inf).*rand(1,model.N)+model.U.inf)']; %check this
% end

while i < max_train_size && falsified==false
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
        x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
        u = [linspace(0,model.T-model.T/model.N,model.N)',((model.U.sup-model.U.inf).*rand(1,model.N)+model.U.inf)']; %check this
    else
        x0=crit_x(1,:)';
        u=crit_u;
    end
    disp(crit_u)

    for j = 1:size(model.spec,1)
    % different types of specifications
        if strcmp(model.spec(j,1).type,'unsafeSet')
            check = any(model.spec(j,1).set.contains(crit_x'));
        elseif strcmp(model.spec(j,1).type,'logic')
            check = ~checkStl(model.spec(j,1).set,crit_x,vpa(linspace(0,model.T,size(crit_x,1)')));
        end
        if check 
            falsified = true;
        end
    end
    i=i+1;
    disp(['iteration completed: ',num2str(i)])
end

%close simulink model
close_system;

train_iter = i;
