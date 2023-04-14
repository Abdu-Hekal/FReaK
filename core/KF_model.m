classdef KF_model
    %Koopman falsification model
    properties
        sim %name of the simulink model. TODO: blackbox function handle
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)
        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)

        T %time horizon for simulation
        dt %time step
        cp %control points for each input signal. needs to be an array of length equal to number of inputs. Needs to be a factor of T/dt
        %default is cp every dt. Note that values other than default are
        %currently only supported for spec of type stl formula.

        soln %internally defined property that stores the solutions for last iteration (do not change)
        specSolns %internally defined property that stores the solutions for each spec for last iteration (do not change)
        cpBool %internal property that is used to set control inputs for pulse inputs (do not change)

        %settings
        maxTrainSize %maximum number of simulations for training before terminating (default: 100)
        trainRand %int, set to 2 to train with random trajectory or 0 to train with previously found crit trajectory or 1 to alternate. (default: 1)
        pulseInput %boolean, set to true if the inputs are pulse inputs, otherwise input is piecewise-constant (default: false)

    end
    methods
        % Constructor
        function model = KF_model(sim)
            model.sim=sim;
            model.maxTrainSize=100;
            model.trainRand=1;
            model.pulseInput = false;
        end

        %simulate function for model, a custom simulate function can be
        %used for a subclass of this class. Ensure that the outputs are consistent
        function [tout, yout, model] = simulate(model, x0, u)
            if isa(model.sim, 'string') || isa(model.sim,"char")
                [tout, yout] = run_simulink(model.sim, model.T, model.dt, x0, u);
            elseif isa(model.sim,'function_handle')
                error('sim as a function handle not yet implemented')
            else
                error('sim not supported')
            end
            model.soln.sims = model.soln.sims+1;
        end

        function [model,trainset]=falsify(model)
            [model,trainset] = coreFalsify(model);
        end
    end
end