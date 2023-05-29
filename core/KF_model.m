classdef KF_model
    %Koopman falsification model
    properties
        model %name of the simulink model. TODO: blackbox function handle
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)
        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)

        T %time horizon for simulation
        dt %time step
        cp %control points for each input signal. needs to be an array of length equal to number of inputs. Needs to be a factor of T/dt
        %default is cp every dt. Note that values other than default are
        %currently only supported for spec of type stl formula.

        %settings
        maxTrainSize %maximum number of simulations for training before terminating (default: 100)
        trainRand %int, set to 3 to train with random trajectory, 2 to train with random neighborhood trajectory, 0 to train with previously found crit trajectory or 1 to alternate between prev and random. (default: 0)
        refineIter %bool, set true to refine with offset. (default: true)
        pulseInput %boolean, set to true if the inputs are pulse inputs, otherwise input is piecewise-constant (default: false)

        %autokoopman settings (struct)
        ak
        %         .obsType: type of observables (default="rff")
        %         .nObs: number of observables (default=100)
        %         .gridSlices: number of slices for grid parameter search (default=5)
        %         .opt: tuner of type "grid", "bopt", or "monte-carlo" (default=grid)
        %         .rank: set of ranks to try of DMD rank parameter (default=[1,200,20]) 

        %solver (optimizer) options (sdpsettings)
        solverOpts

        %internal properties
        soln %internally defined property that stores the solution for last iteration (do not change)
        bestSoln %internally defined property that stores the best solution (do not change)
        specSolns %internally defined property that stores the solutions for each spec for last iteration (do not change)
        cpBool %internal property that is used to set control inputs for pulse inputs (do not change)

    end
    methods
        % Constructor
        function obj = KF_model(model)
            obj.model=model;
            obj.maxTrainSize=100;
            obj.trainRand=0;
            obj.refineIter=true;
            obj.pulseInput = false;

            % autokoopman settings
            obj.ak.obsType="rff";
            obj.ak.nObs=5;
            obj.ak.gridSlices=5;
            obj.ak.opt="grid";
            obj.ak.rank=[0,40,20];

            %default optimizer options
            solver = 'gurobi';  % gurobi, cplex, glpk
            timeLimit = 2000; %2000;
            gapLimit = 100; %0.01;
            gapAbsLimit = inf; %0.1;
            solnLimit = Inf;
            verb = 2;
            obj.solverOpts = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'usex0', 0 ...
                );
            %                 'gurobi.MIPFocus',3,...
            %                 'gurobi.ScaleFlag', 2,...
            %                 'gurobi.BarHomogeneous', 1,...
            %                 'gurobi.CrossoverBasis',-1,...
            %                 'gurobi.DualReductions', 0,...
            %                'cachesolvers',1,...
            %                 'gurobi.InputFile', 'alpha.sol' ...
        end

        %simulate function for model, a custom simulate function can be
        %used for a subclass of this class. Ensure that the outputs are consistent
        function [tout, yout, obj] = simulate(obj, x0, u)
            if isa(obj.model, 'string') || isa(obj.model,"char")
                [tout, yout] = run_simulink(obj.model, obj.T, obj.dt, x0, u);
            elseif isa(obj.model,'function_handle')
                error('sim as a function handle not yet implemented')
            else
                error('sim not supported')
            end
            obj.soln.sims = obj.soln.sims+1;
        end

        function [obj,trainset]=falsify(obj)
            [obj,trainset] = coreFalsify(obj);
        end
    end
end