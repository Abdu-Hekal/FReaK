function [kfModel, A, B,g] = symbolicRFF(kfModel, trainset)

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
    g_ = @(x) sqrt(2/kfModel.ak.nObs)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    path = fileparts(which(mfilename()));
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'autokoopman'));

    g = @(x) autokoopman(x);

end