classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        koopdt %koopman time_step
        solverdt %solver time step for setting up stl, i.e. points where to evaluate stl

        stlList %stl list to falsify
        cpBool %boolean array representing cp points

        %milp sdpvars and constraints
        Falpha %constraints on alpha
        Fstl %constraints for stl
        Freach %constraints states according to reachable sets
        X %states optim var
        Alpha
        Pstl %spdvar storing robustness of stl
        Ostl %spdvar parameters for offset of inequalities in stl

        %optimizer object
        optimizer

        %offset and offset count to modify stl based on prev iterations
        offset
        offsetCount
    end

    properties (Dependent)
        bigM %bigM value for milp
        xlabel %default label name bluSTL
        nx %number of state variables
        L %number of control points
        normz %normalization values based on boundaries of reachable sets
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reachZonos,koopdt,solverdt)
            Sys.reachZonos = reachZonos;
            Sys.koopdt=koopdt;
            Sys.solverdt=solverdt;

            %intitalise offset and count to zero
            Sys.offset=0;
            Sys.offsetCount=-1;
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
            output = {Sys.X,Sys.Alpha,Sys.Pstl};
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
                if Sys.offsetCount>0 %we have an offset
                    param(Sys.offsetCount) = Sys.offset;
                end
                [sol_control, errorflag1,~,~,P] = Sys.optimizer{{param}}; %% call solver
                Sys.X = double(sol_control{1});
                Sys.Alpha = double(sol_control{2});
                Sys.Pstl = double(sol_control{3});
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
        function normz = get.normz(Sys)
            nx=Sys.nx;
            minBound=inf(nx,1);
            maxBound=-inf(nx,1);
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                for k=1:nx
                    supFun=zeros(1,nx);
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

