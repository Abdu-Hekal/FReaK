classdef Koopman_lti

    % Koopman reachability properties
    properties
        reach_zonos = []
        plot = false %plot traces

        xlabel %default label name bluSTL
        nx %number of state variables

        dt %time_step
        L %number of control points
        stl_list %stl_list to falsify

        bigM %
        solver_options %milp solver settings
    end

    methods
        % Constructor
        function Sys = Koopman_lti(reach_zonos,dt)
            Sys.reach_zonos = reach_zonos;
            assert(reach_zonos{1}.isInterval, "initial set must be an interval")
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

            %default bigM
            Sys.bigM=1e6;

        end

        function milp = setup_milp(Sys)
            milp = reach_setup_milp(Sys);
        end
    end
end


