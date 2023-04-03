classdef KF_model
    %Koopman falsification model
    properties
        name %name of the simulink model
        R0 %initial set (CORA class interval)
        U %set of admissible control inputs (class:interval or Zonotope)
        spec %specification defined as an object of the CORA specification class (safe/unsafe sets)

        %bluSTL
        time %time for the dynamics 
        ts %sampling time for controller
        L %no. of (piecewise constant) control inputs

    end
end