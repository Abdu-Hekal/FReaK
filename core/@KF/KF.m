classdef KF
    % KF - class representing a koopman falsification object
    %
    % Syntax:
    %    obj = KF(model)
    %
    % Inputs:
    %    model - name of the simulink model or blackbox function.
    %
    % Outputs:
    %    obj - generated koopman falsification object
    %
    % Example:
    %    obj = KF('Autotrans_shift')
    %    .... %add model properties (see below)
    %    obj.falsify() %falsify model for a given stl
    %
    % See also: KoopMILP

    % Author:      Abdelrahman Hekal
    % Written:      19-November-2023
    % Last update:  ---
    % Last revision:---

    %------------- BEGIN CODE --------------
    properties
        % model: name of the simulink model, blackbox function, OdeFcn or
        % ode (note: ode is a builtin class from matlab 2023b)
        % The blackbox model should be a function handle of the form
        % [tout, yout]=fcn(T,x0,u), where:
        % tout: array of time points for the simulation
        % yout: array of trajectory points corresponding to tout
        % T: time horizon of simulation
        % x0: initial set
        % u: array of inputs, where first column is time points
        model
        % R0: initial set (CORA class interval)
        R0
        % U: set of admissible control inputs (CORA class interval)
        U
        % spec: specification defined as an object of the CORA
        % specification class (safe/unsafe sets/stl)
        spec
        %time horizon for simulation
        T
        % cp: number of control points for each input signal.
        % Needs to be an array of length equal to number of inputs.
        % Needs to be a factor of T/ak.dt. (see below)
        % Default is cp every ak.dt. Note that values other than default
        % are currently only supported for spec of type stl formula.
        cp

        % AutoKoopman settings (struct)
        ak
        %  .dt: koopman time step. default ak.dt=dt.
        %  Change to use coarser koopman step for quicker solution.
        %  Must be a factor of time horizon T.
        %  .obsType: type of observables (default="rff")
        %  .nObs: number of observables (default=100)
        %  .gridSlices: number of slices for grid parameter search (default=5)
        %  .opt: tuner of type "grid", "bopt", or "monte-carlo" (default=grid)
        %  .rank: set of ranks to try of DMD rank parameter (default=[1,200,20])

        % Reachability settings (struct)
        reach
        %   .on: : bool, set to true to use reachability for encoding of MILP.
        % if set to false then direct encoding of the evolution of the
        % koopman linear system as x_{t+1} = A*x_t+B*u_t (default=true).
        % Note that for system's with uncertain initial state, reachability
        % must be used.
        %   .tayOrder: Order of taylor models for reachability (default=6)

        %solver/optimizer (struct)
        solver
        %   .timePoints: solver time points for encoding stl robustness.
        % default solver.timePoints=0:ak.dt:T, i.e. time point every
        % autokoopman step. Change to use coarser solver step when setting
        % up stl constraints for quicker solving time. Use 'auto' to
        % automatically select timePoints using iterative method.
        % Must be no finer than ak.dt
        %   .opts: solver options (see sdpsettings)
        %   .normalize: bool, set to true to normalize optimization objective
        % in milp solver using reachable set bounds, (default=false)
        %   .useOptimizer: bool set to true to use optimizer object. Not
        % using optimizer means stl needs to be setup for milp everytime
        % for offset. setting up optimizer object also takes time.
        % Time trade off? (default=true)

        %SETTINGS
        %runs: number of falsification attempts (default=1)
        runs
        % maxSims: maximum number of simulations for training before
        % terminating, (default=100)
        maxSims
        %timeout: maximum time before algorithm terminates, (default=inf)
        timeout
        % nResets: reset training set after n trajectories (default=5),
        % note that we also reset if milp fails to solve (model is bad)
        nResets
        % trainStrat: int, set to 2 to train with new trajectory, 1 to
        % train with strategy where samples are only addeded to
        %trainset if they improve the robustness,
        % or 0 to train with previously found crit trajectory (default=0)
        trainStrat
        % resetStrat: int, set to 2 to reset to mopso (staliro soar strategy) 
        % trajectory, 1 to reset to perturbed neighborhood trajectory, 
        % or 0 to reset to random trajectory (default=0)
        resetStrat
        % rmRand: bool, set to true to remove first random trajectory when
        % training or false otherwise, (default=true)
        rmRand
        % perturbation percentage increase for generating neighborhood samples
        sampPerturb
        % offsetStrat: int, set 1 to refine with offset, 0 for no offset,
        % -1 to offset next iteration (after retraining koopman model),
        % (default=-1).
        offsetStrat
        % dt: time step. This time step is used to interpolate inputs to
        % system . Use a small dt for accurate results, (default:0.01)
        % Must be factor of time horizon T.
        dt
        % Interpolation types for input & trajectory. See "interp1" for supported types
        % inputInterpolation: interpolate input between control points.
        % (default='previous'). Note that control points may decrease
        % with coarser time step for koopman.
        inputInterpolation
        % trajInterpolation: interpolate output trajectory for learning
        % autokoopman model, (default='linear')
        trajInterpolation
        % pulseInput: boolean, set to true if the inputs are pulse inputs,
        % otherwise false (default=false)
        pulseInput
        % verb: int, set 3 for printing all info whilst falsifying, 
        % 2 for only best soln info whilst falsifying, 1 for
        % print at end of falsifying, 0 for no print (default=0)
        verb

        %internal properties (DO NOT CHANGE)
        cpBool %used to set control inputs for pulse inputs
        inputsInterval %1d interval of all applicable inputs (x0 and/or u)


    end
    methods
        % Constructor
        function obj = KF(model)
            obj.model=model;
            obj.runs=1;
            obj.dt=0.01;
            obj.maxSims=5000;
            obj.timeout=inf;
            obj.nResets=5; %5
            obj.trainStrat=0;
            obj.resetStrat=0;
            obj.rmRand=true; %true
            obj.sampPerturb=0.05;
            obj.offsetStrat=-1; %-1
            obj.inputInterpolation='previous';
            obj.trajInterpolation='linear';
            obj.pulseInput = false;
            obj.verb=0;

            % autokoopman settings
            obj.ak.obsType="rff";
            obj.ak.nObs=20;
            obj.ak.gridSlices=5;
            obj.ak.opt="grid"; %grid
            obj.ak.rank=[1,20,4];

            %reachability settings
            obj.reach.on=true; %true
            obj.reach.tayOrder=6;

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
                'gurobi.MIPFocus',3,...
                'gurobi.NumericFocus',2,...
                'gurobi.BarHomogeneous', 1,...
                'usex0', 0 ...
                );
            obj.solver.normalize=false;
            obj.solver.useOptimizer=true;
        end
    end
end