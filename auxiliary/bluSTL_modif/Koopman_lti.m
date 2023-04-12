classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        dt %time_step

        xlabel %default label name bluSTL
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
        nx %number of state variables
        L %number of control points
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reachZonos,dt)
            Sys.reachZonos = reachZonos;
            Sys.dt=dt;

            % default label names
            Sys.xlabel = cell(1,Sys.nx);
            for iX = 1:Sys.nx
                Sys.xlabel{iX} = ['x' num2str(iX)];
            end

            %default solver options:
            solver = 'gurobi';  % gurobi, cplex, glpk
            timeLimit = 2000;
            gapLimit = 0.01;
            gapAbsLimit = 0.1;
            solnLimit = Inf;
            verb = 0;
            Sys.solver_options = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'cachesolvers',1,...
                'gurobi.BarHomogeneous', 1,...
                'gurobi.ScaleFlag', 2, ...
                'usex0',1, ...
                'gurobi.DualReductions', 0);
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

        function value = get.nx(Sys)
            value=size(Sys.reachZonos{1}.center,1);
        end

        function value = get.L(Sys)
            value=size(Sys.reachZonos,1)-1;
        end

        function value = get.bigM(Sys)
            %find suitable bigM based on zonotope boundaries
            value=0;
            for i=1:length(Sys.reachZonos)
                zono=Sys.reachZonos{i};
                if norm(zono,inf) > value
                    order = ceil(log10(norm(zono)));
                    value = 10^order;
                end
            end
        end

        function diagnostics = optimize(Sys)

            constraints=[Sys.Falpha, Sys.Fstl, Sys.Freach];
            objective=Sys.Pstl; %objective is to miniimize robustness
            options=Sys.solver_options;
            %% call solverarch
            diagnostics = optimize(constraints,objective,options);
        end
    end
end

