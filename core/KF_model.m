classdef KF_model
    %Koopman falsification model
    properties
        model %name of the simulink model. TODO: blackbox function handle
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)
        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)

        T %time horizon for simulation
        dt %time step. This time step is used to interpolate inputs to system and output trajecotries. Use a small dt for accurate results, (default:0.01)
        cp %max control points for each input signal. needs to be an array of length equal to number of inputs. Needs to be a factor of T/dt & T/kdt
        %default is cp every dt. Note that values other than default are
        %currently only supported for spec of type stl formula.

        %settings
        maxSims %maximum number of simulations for training before terminating (default: 100)
        nResets %reset training set after n trajectories (default: 20), note that we also reset if milp fails to solve (model is bad)
        trainRand %int, set to 3 to train with random trajectory, 2 to train with random neighborhood trajectory, 0 to train with previously found crit trajectory or 1 to alternate between prev and random. (default: 0)
        offsetStrat %int, set 1 to refine with offset, 0 for no offset, -1 to offset next iteration (after retraining koopman model). (default: 1). Note that offset strategy 1 is most stable and strategy -1 can lead to problems when warmstart (usex0) is used.
        useOptimizer %bool set to true to use optimizer object. Not using optimizer means stl needs to be setup for milp everytime for offset. setting up optimizer object also takes time. Time trade off? Note that using optimzier is most stable and not using can lead to problems when warmstart (usex0) is used.
        reach %use reachability for encoding of MILP (default:true)
        % interpolation types for input & trajectory. See "interp1" for supported types
        inputInterpolation % interpolate input between control points. default 'previous'. Note that control points may decrease with coarser koopman
        trajInterpolation %interpolate output trajectory for learning autokoopman model, default 'linear'
        pulseInput %boolean, set to true if the inputs are pulse inputs, otherwise input is piecewise-constant (default: false)

        %autokoopman settings (struct)
        ak
        %          .dt: koopman time step. default ak.dt=dt. Change to use coarser koopman step for quicker solution
        %         .obsType: type of observables (default="rff")
        %         .nObs: number of observables (default=100)
        %         .gridSlices: number of slices for grid parameter search (default=5)
        %         .opt: tuner of type "grid", "bopt", or "monte-carlo" (default=grid)
        %         .rank: set of ranks to try of DMD rank parameter (default=[1,200,20])

        %solver/optimizer (struct)
        solver
        %         .dt: solver time step. default solver.dt=ak.dt Change to use coarser solver step when setting up stl constraints for quicker solution
        %         .opts: solver options (sdpsettings)

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
            obj.dt=0.01; 
            obj.maxSims=5000;
            obj.nResets=5;
            obj.trainRand=0;
            obj.offsetStrat=-1;
            obj.useOptimizer=false;
            obj.reach=true;
            obj.inputInterpolation='previous';
            obj.trajInterpolation='linear';
            obj.pulseInput = false;
            obj.soln=struct;

            % autokoopman settings
            obj.ak.obsType="rff";
            obj.ak.nObs=20;
            obj.ak.gridSlices=5;
            obj.ak.opt="grid"; %grid
            obj.ak.rank=[1,20,4];

            %default optimizer options
            solver = 'gurobi';  % gurobi, cplex, glpk
            timeLimit = 120; %2000;
            gapLimit = 10e-4; %0.1;
            gapAbsLimit = 10e-10; %0.1;
            solnLimit = Inf;
            verb = 0;
            obj.solver.opts = sdpsettings('verbose', verb,'solver', solver, ...
                'gurobi.TimeLimit', timeLimit, ...
                'gurobi.MIPGap', gapLimit, ...
                'gurobi.MIPGapAbs', gapAbsLimit, ...
                'gurobi.SolutionLimit', solnLimit,...
                'gurobi.Method',3,...
                'cachesolvers',1,...
                'usex0', 0 ...
                );

%                 'gurobi.MIPFocus',3,...
%                 'gurobi.NumericFocus',2,...
%                 'gurobi.BarHomogeneous', 1,...
%                 'cachesolvers',1,...

            %                'gurobi.NumericFocus',3,...
            %                 'gurobi.MIPFocus',3,...
            %                 'gurobi.ScaleFlag', 2,...
            %                 'gurobi.BarHomogeneous', 1,...
            %                 'gurobi.CrossoverBasis',-1,...
            %                 'gurobi.DualReductions', 0,...
            %                'cachesolvers',1,...
            %                 'gurobi.InputFile', 'alpha.sol' ...


            %create empty struct to store prev soln
            obj.soln=struct;
            obj.soln.koopTime=0; obj.soln.milpSetupTime=0; obj.soln.milpSolvTime=0; obj.soln.simTime=0;
            obj.soln.sims=0;
            %create empty struct to store best soln
            obj.bestSoln=struct; obj.bestSoln.rob=inf; obj.bestSoln.timeRob=inf;
            %create empty dict to store prev soln for each spec
            obj.specSolns = dictionary(obj.spec,struct);
        end

        %simulate function for model, a custom simulate function can be
        %used for a subclass of this class. Ensure that the outputs are consistent
        function [tout, yout, obj] = simulate(obj, x0, u)
            tic
            if isa(obj.model, 'string') || isa(obj.model,"char")
                [tout, yout] = run_simulink(obj.model, obj.T, x0, u);
            elseif isa(obj.model,'function_handle')
                %function handle must have 3 inputs T,x0,u
                [tout, yout] = obj.model(obj.T, x0, u);
            else
                error('model not supported')
            end
            obj.soln.sims = obj.soln.sims+1;
            obj.soln.simTime = obj.soln.simTime+toc;
        end

        function [tout, yout, x0, u] = randSimulation(obj)
            [x0,u] = getRandomSampleXU(obj);
            [tout, yout,~] = simulate(obj, x0, u);
            end

        function [obj,trainset]=falsify(obj)
            [obj,trainset] = coreFalsify(obj);
        end
    end
end