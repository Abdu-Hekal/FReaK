classdef KF_model
    %Koopman falsification model
    properties
        sim %name of the simulink model. TODO: blackbox function handle
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)

        T %time horizon for simulation
        dt %time step
        %TODO: varying time control points.
        cp %control points for each input signal. needs to be an array of length equal to number of inputs. Needs to be a factor of T/dt
        %default is cp every dt. Note that values other than default are
        %currently only supported for spec of type stl formula.

        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)
        spec_soln %internally defined property that stores the solutions for each spec (do not change)

        pulse_input %boolean, set to true if the inputs are pulse inputs, otherwise input is piecewise-constant (default: false)
        cp_bool %internal property that is used to set control inputs for pulse inputs (do not change)

    end
    methods
        % Constructor
        function model = KF_model(sim)
            model.sim=sim;
            model.pulse_input = false;
        end

        %simulate function for model, a custom simulate function can be
        %used for a subclass of this class. Ensure that the outputs are consistent 
        function [tout, yout] = simulate(model, x0, u)
            if isa(model.sim, 'string') || isa(model.sim,"char")
                [tout, yout] = run_simulink(model.sim, model.T, model.dt, x0, u);
            elseif isa(model.sim,'function_handle')
                error('sim as a function handle not yet implemented')
            else
                error('sim not supported')
            end
        end
    end
end