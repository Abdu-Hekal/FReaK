function [sys, X, U, t, x, crit_x] = modelvanderpol(X, U, t, T)

    %generate training data
    x = 0.3*rand()+1.25;
    y = 0.1*rand()+2.25;
    init_set = [x,y];
    [t{end+1}, x] = run_vanderpol(init_set, T);
    X{end+1} = x';
    U{end+1} = zeros(size(x')-1);
    dt = t{1}(2) - t{1}(1);
    
    % identify Koopman model
    numFeat = 100;                      % number observables
    rank = 5;                           % rank for DMD
    l = 0.002;                          % lengthscale for kernel
    
    [cost,g,A,B,C,c] = costFunFourier([numFeat;rank;l;1],X,U,dt,false,false);  
    
    % save observables in function
    path = fileparts(which(mfilename()));
    xSym = sym('x',[size(X{1},1),1]);
    matlabFunction(g,'Vars',{xSym},'File',fullfile(path,'functionVanderPol'));
    
    % visualize the predictions for the identfied Koopman model
    sys = linearSysDT(linearSys(A,B,c,C),dt);
    
    % reachability analysis (reachable set for all possible inputs)
    params.R0 = zonotope(functionVanderPol(interval([1.25;2.25],[1.55;2.35])));
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

