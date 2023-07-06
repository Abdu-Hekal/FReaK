classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        U %set of admissible control inputs (class:interval or Zonotope)
        koopdt %koopman time_step
        solverdt %solver time step for setting up stl, i.e. points where to evaluate stl

        stl %stl to falsify
        cpBool %boolean array representing cp points

        %milp sdpvars and constraints
        Falpha %constraints on alpha
        Fstl %constraints for stl
        Freach %constraints states according to reachable sets
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
        xlabel %variables label name from bluSTL
        ulabel %input label names
        nx %number of state variables
        L %number of control points
        normz %normalization values based on boundaries of reachable sets
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reachZonos,U,koopdt,solverdt)
            Sys.reachZonos = reachZonos;
            Sys.U = U;
            Sys.koopdt=koopdt;
            Sys.solverdt=solverdt;

            %intitalise offset map to empty 
            Sys.offsetMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
        end

        function Sys = setupAlpha(Sys)
            Sys = KoopmanSetupAlpha(Sys);
        end

        function Sys = setupStl(Sys,hardcoded)
            Sys = koopmanSetupStl(Sys,hardcoded);
        end

        function Sys = setupReach(Sys)
            Sys = KoopmanSetupReach(Sys);
        end

        function Sys = setupOptimizer(Sys,options)
            constraints=[Sys.Falpha, Sys.Fstl, Sys.Freach];
            objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
            output = {Sys.x,Sys.alpha,Sys.u,Sys.Pstl};
            % setup optimizer
            Sys.optimizer = optimizer(constraints,objective,options,Sys.Ostl, output);
        end

        function Sys = optimize(Sys,options)
            if isempty(Sys.optimizer) %no optimizer object, optimize directly
                constraints=[Sys.Falpha, Sys.Fstl, Sys.Freach];
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
                assign(Sys.alpha,double(sol_control{2}));
                assign(Sys.u,double(sol_control{3}));
                assign(Sys.Pstl,double(sol_control{4}));
                if errorflag1 ~= 0
                   disp(['Yalmip error: ' yalmiperror(errorflag1)]); % some other error
                end
            end
        end

        %getters for dependent properties
        function bigM = get.bigM(Sys)
            %find suitable bigM based on zonotope boundaries
            bigM=0;
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                if norm(zono,inf) > bigM
                    order = ceil(log10(norm(zono)));
                    bigM = 10^order;
                end
            end
        end
        function nx = get.nx(Sys)
            nx=size(Sys.reachZonos{1}.center,1);
        end
        function L=get.L(Sys)
            L=size(Sys.reachZonos,1)-1;

        end
        function xlabel=get.xlabel(Sys)
            % default label names
            xlabel = cell(1,Sys.nx);
            for iX = 1:Sys.nx
                xlabel{iX} = ['x' num2str(iX)];
            end
        end
        function ulabel=get.ulabel(Sys)
            % default label names
            ulabel = cell(1,size(Sys.U,1));
            for iU = 1:size(Sys.U,1)
                ulabel{iU} = ['u' num2str(iU)];
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

