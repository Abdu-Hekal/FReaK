function [falsified, trainset, crit_x, train_iter] = coreFalsify(model, max_train_size)

%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

%generate random initial set and input
x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
u = [linspace(0,model.T-model.T/model.N,model.N)',((model.U.sup-model.U.inf).*rand(1,model.N)+model.U.inf)']; %check this

trainset.X = {}; trainset.XU={}; trainset.t = {}; %empty cells to store states, inputs and times for training trajectories

falsified = false;
i = 0;
while i < max_train_size && falsified==false
%     disp("input")
%     disp(u)
    [trainset, crit_x, crit_u] = symbolicRFF(model, trainset, x0, u);
    %retrain with initial set as the critical set found in prev iteration
    x0=crit_x(1,:)';
    u=crit_u;
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
