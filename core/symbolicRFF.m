function [kfModel, trainset] = symbolicRFF(kfModel, trainset, x0, u)

    if size(x0,2) > 1
        x=x0; %x0 passed is infact the full critical trajectory
        trainset.t{end+1} = trainset.t{end}; %time step is same anyway
    else
        [trainset.t{end+1}, x, kfModel] = simulate(kfModel, x0, u);
    end
    trainset.X{end+1} = x';
    %append zeros at end to account for last time point (which has no
    %inputs), but length must be consistant with trajectory states
    trainset.XU{end+1} = [u(:,2:end)', zeros(size(u,2)-1,1)];

    % AutoKoopman settings and run
    if ~isempty(kfModel.U) %check if kfModel has inputs
        inputs_list = trainset.XU;
    else
        inputs_list = string(missing);
    end 
    tic
    pyrunfile("run_autokoopman.py",'koopman_model',times=trainset.t,trajectories=trainset.X, param_dict=kfModel.ak,inputs_list=inputs_list);
    kfModel.soln.koopTime= kfModel.soln.koopTime+toc;
    load("autokoopman_model.mat", "A","B","u","w")

    % create observables function
    n = size(kfModel.R0,1); %number of variables
    g_ = @(x) sqrt(2/20)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    path = fileparts(which(mfilename()));
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'autokoopman'));

    g = @(x) autokoopman(x);
    [crit_x0,crit_u, kfModel] = falsifyFixedModel(A,B,g,kfModel);
    all_steps = kfModel.T/kfModel.dt;
    crit_u = [crit_u';zeros(size(crit_u,1),all_steps-size(crit_u,2))'];
    kfModel.soln.u = [linspace(0,kfModel.T-kfModel.dt,all_steps)',crit_u];

    % run most critical input on the real system
    [t, x, kfModel] = simulate(kfModel, crit_x0, kfModel.soln.u);
    kfModel.soln.t=t;
    kfModel.soln.x=x;
end