function R = reachabilityAnalysis(sys,g,x0,u,tFinal,V,W)
% compute the reachable set for the Koopman operator linearized system

    % modified linear system
    sys = linearSysDT(sys.A,[sys.B eye(size(sys.A,1))],sys.c,sys.C,sys.dt);

    % reachability parameter
    params.R0 = zonotope(g(x0));
    params.U = cartProd(zeros(size(u,1),1),zonotope(W));
    params.tFinal = tFinal;
    params.u = [u;zeros(size(sys.A,1),size(u,2))];

    % reachabiltiy settings
    options.zonotopeOrder = 20;

    % reachabiltiy analysis
    R = reach(sys,params,options);
    R = R + V;
end