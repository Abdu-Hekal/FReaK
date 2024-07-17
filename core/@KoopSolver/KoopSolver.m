classdef KoopSolver
% KoopSolver - class representing a koopman solver object for setting up and
%   solving optimization problem which minimizes robustness of koopman model w.r.t stl.
%
% Syntax:
%    obj = KoopSolver(T,koopdt,solverTimePoints,X0,U)
%
% Inputs:
%    T - time horizon for koopman model.
%    koopdt - time step for koopman model
%    solverTimePoints - time points for stl analysis (multiple of koopdt)
%    X0 - initial set (class:interval or Zonotope)
%    U - set of admissible control inputs (class:interval or Zonotope)
%
% Outputs:
%    obj - generated koopman falsification object
%
%
% See also: KF

% Author:      Abdelrahman Hekal
% Written:      19-November-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
properties
    reachZonos = [] %reachable zonotopes, used if reachable encoding.
    % A and B matrices in linear evolution of system, Ax+Bu, only used
    %   if direct encoding.
    A
    B
    g % observables function,
    X0 %initial set (not that it should be a point if reachability is not used)
    U %set of admissible control inputs (class:interval or Zonotope)
    koopdt %koopman time_step
    solverTimePoints %solver time points for setting up stl, i.e. points where to evaluate stl
    nx %number of state variables
    nu %number of inputs
    nObs %number of observables
    L %number of control points
    normalize %bool set to true to normalize optimization objective using reachable set bounds

    stl %stl to falsify
    cpBool %boolean array representing cp points

    %sdpvars and constraints
    Finit %constraints on alpha or inputs
    Fstl %constraints for stl
    Fdyn %constraints on states
    x %states optim var
    alpha %alpha optim var
    u %inputs optim var
    Pstl %spdvar storing robustness of stl
    Ostl %spdvar parameters for offset of inequalities in stl (when using optimizer object)
    Wstl %spdvar parameters for weights in weighted stl (when using optimizer object)

    %optimizer object
    optimizer

    %offset and offset count to modify stl based on prev iterations
    offsetMap
    %weights for weighted stl (recusrive manipulation of weights to falsify stl)
    weights
    % maximum weight to ensure that problem does not become too large (default 10e9)
    maxWeight
end

properties (Dependent)
    bigM %bigM value for milp
    normz %normalization values based on boundaries of reachable sets
end

methods
    % Constructor
    function Sys = KoopSolver(T,koopdt,solverTimePoints,X0,U)
        Sys.L=ceil(T/koopdt);
        Sys.koopdt=koopdt;
        Sys.solverTimePoints=solverTimePoints;

        Sys.X0 = X0;
        Sys.U = U;

        Sys.nx=size(Sys.X0,1);
        Sys.nu=size(Sys.U,1);

        %intitalise offset map to empty
        Sys.offsetMap = dictionary();
        %initialise weights to empty array 
        Sys.weights = [];
        % default not to use normalization
        Sys.normalize=false;
        %setup (unconstrained) objective fcn
        Sys.Pstl=sdpvar(1,1);
        %default max weight
        Sys.maxWeight=10e9;
    end

    %getters for dependent properties
    function bigM = get.bigM(Sys)
        if ~isempty(Sys.reachZonos)
            %find suitable bigM based on zonotope boundaries
            bigM=0;
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                if norm(zono,inf) > bigM
                    order = ceil(log10(norm(zono)));
                    bigM = 10^order;
                end
            end
            %if input is not empty
            if ~isempty(Sys.U)
                %make sure bigM is bigger also than inputs
                inpmax = 10^(ceil(log10(max(Sys.U.sup))));
                bigM=max(bigM,inpmax);
            end
        else
            bigM = 10e6;
        end
        %sanity check, make sure bigM is not too large
        bigM = min(bigM,10e6);
    end
    function normz = get.normz(Sys)
        minBound=inf(Sys.nx,1);
        maxBound=-inf(Sys.nx,1);
        for i=1:length(Sys.reachZonos)
            zono=Sys.reachZonos{i};
            for k=1:Sys.nx
                supFun=zeros(1,Sys.nx);
                supFun(k)=-1;
                minBound(k) = min(minBound(k),-supportFunc(zono,supFun));
                supFun(k)=1;
                maxBound(k) = max(maxBound(k),supportFunc(zono,supFun));
            end
        end
        normz = maxBound-minBound;
    end
end
end

