classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        dt %time_step

        stlList %stl list to falsify
        cpBool %boolean array representing cp points

        solver_options %milp solver settings

        %milp sdpvars and constraints
        Falpha %constraints on alpha
        Fstl %constraints for stl
        Freach %constraints states according to reachable sets
        X %states optim var
        Alpha
        Pstl %robustness of stl

        %optimizer object
        milp
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
        function Sys = Koopman_lti(reachZonos,dt)
            Sys.reachZonos = reachZonos;
            Sys.dt=dt;

            %default solver options:
            solver = 'gurobi';  % gurobi, cplex, glpk
            timeLimit = 2000; %2000;
            gapLimit = 0.00001; %0.01;
            gapAbsLimit = inf; %0.1;
            solnLimit = Inf;
            verb = 1;
            Sys.solver_options = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'usex0', 1 ...
);
%                 'gurobi.ScaleFlag', 2,...
%                 'gurobi.BarHomogeneous', 1,...
%                 'gurobi.CrossoverBasis',-1,...
%                 'gurobi.DualReductions', 0,...
%                'cachesolvers',1,...
%                 'gurobi.InputFile', 'alpha.sol' ... 

        end

        function Sys = setupAlpha(Sys)
            Sys = KoopmanSetupAlpha(Sys);
        end

        function Sys = setupStl(Sys)
            Sys = koopmanSetupStl(Sys);
        end

        function Sys = setupReach(Sys)
            Sys = KoopmanSetupReach(Sys);
        end

        function diagnostics = optimize(Sys)

            constraints=[Sys.Falpha, Sys.Fstl, Sys.Freach];
            objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
            options=Sys.solver_options;
            %% call solverarch
            diagnostics = optimize(constraints,objective,options);
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

