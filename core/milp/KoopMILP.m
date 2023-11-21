classdef KoopMILP
% KoopMILP - class representing a koopman milp object for setting up and
%   solving MILP which minimizes robustness of koopman model w.r.t stl.
%
% Syntax:
%    obj = KoopMILP(T,koopdt,solverdt,X0,U)
%
% Inputs:
%    T - time horizon for koopman model.
%    koopdt - time step for koopman model
%    solverdt - time step for stl analysis (multiple of koopdt)
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
    solverdt %solver time step for setting up stl, i.e. points where to evaluate stl
    nx %number of state variables
    nu %number of inputs
    nObs %number of observables
    L %number of control points
    normalize %bool set to true to normalize optimization objective using reachable set bounds

    stl %stl to falsify
    cpBool %boolean array representing cp points

    %milp sdpvars and constraints
    Finit %constraints on alpha or inputs
    Fstl %constraints for stl
    Fdyn %constraints on states
    x %states optim var
    alpha %alpha optim var
    u %inputs optim var
    Pstl %spdvar storing robustness of stl
    Ostl %spdvar parameters for offset of inequalities in stl

    %optimizer object
    optimizer

    %offset and offset count to modify stl based on prev iterations
    offsetMap
end

properties (Dependent)
    bigM %bigM value for milp
    normz %normalization values based on boundaries of reachable sets
end

methods
    % Constructor
    function Sys = KoopMILP(T,koopdt,solverdt,X0,U)
        Sys.L=ceil(T/koopdt);
        Sys.koopdt=koopdt;
        Sys.solverdt=solverdt;

        Sys.X0 = X0;
        Sys.U = U;

        Sys.nx=size(X0,1);
        Sys.nu=size(U,1);

        %intitalise offset map to empty
        Sys.offsetMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
        % default not to use normalization
        Sys.normalize=false;
    end

    function Sys = setupAlpha(Sys)
        Sys = koopSetupAlpha(Sys);
    end

    function Sys = setupStl(Sys,hardcoded)
        Sys = koopSetupStl(Sys,hardcoded);
    end

    function Sys = setupReach(Sys)
        Sys = koopSetupReach(Sys);
    end

    function Sys = setupInit(Sys)
        Sys = koopSetupInit(Sys);
    end

    function Sys = setupDynamics(Sys)
        Sys = koopSetupDynamics(Sys);
    end

    function Sys = setupOptimizer(Sys,options)
        constraints=[Sys.Finit, Sys.Fstl, Sys.Fdyn];
        objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
        if ~isempty(Sys.reachZonos)
            if ~isempty(Sys.u)
                output = {Sys.x,Sys.Pstl,Sys.alpha,Sys.u};
            else
                output = {Sys.x,Sys.Pstl,Sys.alpha};
            end
        else
            output = {Sys.x,Sys.Pstl,Sys.u};
        end
        % setup optimizer
        Sys.optimizer = optimizer(constraints,objective,options,Sys.Ostl, output);
    end

    function Sys = optimize(Sys,options)
        if isempty(Sys.optimizer) %no optimizer object, optimize directly
            constraints=[Sys.Finit, Sys.Fstl, Sys.Fdyn];
            objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
            %% call solverarch
            optimize(constraints,objective,options);
        else
            param = zeros(1,length(Sys.Ostl));
            if Sys.offsetMap.Count > 0 %we have an offset
                keys = Sys.offsetMap.keys;
                for ii=1:Sys.offsetMap.Count
                    key = keys{ii};
                    param(key) = Sys.offsetMap(key);
                end
            end
            [sol_control, errorflag1,~,~,P] = Sys.optimizer{{param}}; %% call solver
            assign(Sys.x,double(sol_control{1}));
            assign(Sys.Pstl,double(sol_control{2}));
            if ~isempty(Sys.reachZonos)
                assign(Sys.alpha,double(sol_control{3}));
                if ~isempty(Sys.u)
                    assign(Sys.u,double(sol_control{4}));
                end
            else
                assign(Sys.u,double(sol_control{3}));
            end
        end
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

