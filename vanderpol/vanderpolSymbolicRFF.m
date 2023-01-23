function [sys, X, t, x, crit_x] = vanderpolSymbolicRFF(init_set, X, t, T)
    

    [t{end+1}, x] = run_vanderpol(init_set, T);
    X{end+1} = x';
    dt = t{1}(2) - t{1}(1);

    % AutoKoopman settings and run
    param_dict = struct("samp_period", dt, "obs_type", "rff", "n_obs",100, "grid_param_slices", 5);
    pyrunfile("run_autokoopman.py",'koopman_model',times=t,trajectories=X, param_dict=param_dict,inputs_list=string(missing));
    load("autokoopman_model.mat")

    % create observables function
    g_ = @(x) sqrt(2/100)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    n = 2;
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    path = fileparts(which(mfilename()));
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'autokoopmanVanderPol'));
    
    % visualize the predictions for the identfied Koopman model
    c = zeros(102,1);
    B = zeros(102,1);
    C = [eye(n), zeros(n,100)];
    sys = linearSysDT(A,B,c,C,dt); 

    % reachability analysis (reachable set for all possible inputs)
    params.R0 = zonotope(autokoopmanVanderPol(interval([1.25;2.25],[1.55;2.35])));
    params.U = zonotope(0);
    params.tFinal = 7;  
    options.zonotopeOrder = 10000;
    
    R = reach(sys,params,options);
    
    % determine most critical input from reachable set
    ind = []; val = -inf;
    
    for i = 1:length(R.timePoint.set)
        val_ = supportFunc(R.timePoint.set{i},[1 0]);
        if val_ > val
            val = val_; ind = i;
        end
    end
    
    [~,~,alpha] = supportFunc(R.timePoint.set{ind},[1 0]);
    critical_init_set = center(params.R0) + generators(params.R0) * alpha;
    
    % run most critical input on the real system
    [~, crit_x] = run_vanderpol(critical_init_set, T);

end