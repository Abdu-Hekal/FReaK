classdef KF_model
    %Koopman falsification model
    properties
        sim %name of the simulink model. TODO: blackbox function handle
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)

        T %time horizon for simulation
        dt %time step
        cp %control points for each input signal. needs to be an array of length equal to number of inputs. Needs to be a factor of T/dt
        %default is cp every dt. Note that values other than default are
        %currently only supported for spec of type stl formula.

        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)
        spec_soln %internally defined property that stores the solutions for each spec (do not change)

        pulse_input %boolean, set to true if the inputs are pulse inputs, otherwise input is piecewise-constant
        cp_bool %internal property that is used to set control inputs for pulse inputs (do not change)

    end
    methods
        % Constructor
        function model = KF_model(sim)
            model.sim=sim;
            model.pulse_input = true;
        end
    end
end