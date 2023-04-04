function [trainset, crit_x, crit_u] = symbolicRFF(model, trainset, x0, u)

    [trainset.t{end+1}, x] = run_simulation(model.sim, model.T, x0, u);
    trainset.X{end+1} = x';
    dt = trainset.t{1}(2) - trainset.t{1}(1);
    %append zeros at end to account for last time point (which has no
    %inputs), but length must be consistant with trajectory states
    trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];

    % AutoKoopman settings and run
    param_dict = struct("samp_period", dt, "obs_type", "rff", "n_obs",100, "grid_param_slices", 5);
    if ~isempty(find_system(model.sim,'BlockType','Inport')) %check if simulink model has inputs
        inputs_list = trainset.XU;
    else
        inputs_list = string(missing);
    end    
    pyrunfile("run_autokoopman.py",'koopman_model',times=trainset.t,trajectories=trainset.X, param_dict=param_dict,inputs_list=inputs_list);
    load("autokoopman_model.mat", "A","B","u","w")

    % create observables function
    n = size(model.R0,1); %number of variables
    g_ = @(x) sqrt(2/100)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    path = fileparts(which(mfilename()));
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'autokoopman'));

    g = @(x) autokoopman(x);
    [crit_x0,crit_u] = falsifyFixedModel(A,B,g,dt,model.spec,model.R0,model.U,model.T,model.cp_bool);
    all_steps = model.T/model.dt;
    crit_u = [crit_u';zeros(size(crit_u,1),all_steps-size(crit_u,2))'];
    crit_u = [linspace(0,model.T-model.dt,all_steps)',crit_u];

    % run most critical input on the real system
    [~, crit_x] = run_simulation(model.sim, model.T, crit_x0, crit_u);

end