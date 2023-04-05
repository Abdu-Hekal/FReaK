classdef Koopman_lti

    % Koopman reachability properties
    properties
        reach_zonos = [] %reachable zonotopes
        dt %time_step
        cp_bool %boolean array representing cp points

        xlabel %default label name bluSTL
        nx %number of state variables
        L %number of control points
        stl_list %stl_list to falsify

        bigM %
        solver_options %milp solver settings


        Fstl
        Falpha
        output
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reach_zonos,dt)
            Sys.reach_zonos = reach_zonos;
            Sys.dt=dt;
            Sys.nx=size(reach_zonos{1}.center,1);
            Sys.L=size(reach_zonos,1)-1;

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
                'gurobi.DualReductions', 0);

            %find suitable bigM based on zonotope boundaries
            Sys.bigM=0;
            for i=1:length(reach_zonos)
                zono=reach_zonos{i};
                if norm(zono,inf) > Sys.bigM
                    order = ceil(log10(norm(zono)));
                    Sys.bigM = 10^order;
                end
            end
        end

        function Sys = setup_milp(Sys)
            Sys = koopman_setup_milp(Sys);
        end
        
        function milp = reach_milp(Sys)
            milp = koopman_reach_milp(Sys);
        end
    end
end


