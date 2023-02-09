function [sys, X, t, x, crit_x] = vanderpolSymbolicRFF(init_set, X, t, T)
    

    [t{end+1}, x] = run_vanderpol(init_set, T);
    X{end+1} = x';
    dt = t{1}(2) - t{1}(1);

    % AutoKoopman settings and run
    param_dict = struct("samp_period", dt, "obs_type", "rff", "n_obs",100, "grid_param_slices", 5);
    pyrunfile("run_autokoopman.py",'koopman_model',times=t,trajectories=X, param_dict=param_dict,inputs_list=string(missing));
    load("autokoopman_model.mat", "A","B","u","w")

    % create observables function
    n = 2; %number of variables
    g_ = @(x) sqrt(2/100)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    path = fileparts(which(mfilename()));
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'autokoopmanVanderPol'));
    
    % visualize the predictions for the identfied Koopman model
    if isempty(B)
        B = zeros(size(A,1),1);
    end
    c = zeros(size(A,1),1);
    C = [eye(n), zeros(n,size(A,1)-n)];
    sys = linearSysDT(A,B,c,C,dt); 

    spec = specification(halfspace([1 0],2.095),'unsafeSet');
    g = @(x) autokoopmanVanderPol(x);
    R0 = interval([1.25;2.25],[1.55;2.35]);
    U = interval(0,0);
    tFinal = 7;
    [x0,u] = falsifyFixedModel(A,B,g,dt,spec,R0,U,tFinal);
    
    % run most critical input on the real system
    [~, crit_x] = run_vanderpol(x0, T);

end