classdef Koopman_lti
    properties
        reachZonos = [] %reachable zonotopes
        dt %time_step

        xlabel %default label name bluSTL
        nx %number of state variables
        L %number of control points
        stlList %stl list to falsify
        cpBool %boolean array representing cp points

        bigM %bigM value for milp
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

    methods
        % Constructor
        function Sys = Koopman_lti(reachZonos,dt)
            Sys.reachZonos = reachZonos;
            Sys.dt=dt;
            Sys.nx=size(reachZonos{1}.center,1);
            Sys.L=size(reachZonos,1)-1;

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
            verb = 1;
            Sys.solver_options = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'cachesolvers',1,...
                'gurobi.BarHomogeneous', 1,...
                'gurobi.ScaleFlag', 2, ...
                'usex0',1, ...
                'gurobi.DualReductions', 0, ...
                'gurobi.TuneTrials',0,... %default
                'gurobi.PreSOS2BigM',-1,... %default
                'gurobi.CrossoverBasis',-1); %default
%                 'gurobi.NumericFocus',1,...
            Sys.solver_options = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'cachesolvers',1,...
                'usex0',1);

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

        %reach zonos setter, which also sets bigM value
        function Sys=set.reachZonos(Sys,reachZonos)
            Sys.reachZonos=reachZonos;
            %find suitable bigM based on zonotope boundaries
            Sys.bigM=0;
            for i=1:length(reachZonos)
                zono=reachZonos{i};
                if norm(zono,inf) > Sys.bigM
                    order = ceil(log10(norm(zono)));
                    Sys.bigM = 10^order;
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

