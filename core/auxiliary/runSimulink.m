function [tout, yout,xFinal] = runSimulink(model_name, T, x0, u)
    % runSimulink - Run a simulation of a Simulink model
    %
    % Syntax:
    %    [tout, yout] = runSimulink(model_name, T, x0, u)
    %
    % Description:
    %    This function runs a simulation of a Simulink model, given the model
    %    name, simulation time T, initial state x0, and input signal u. It
    %    extracts the initial state information and updates the simulation
    %    parameters accordingly, then performs the simulation and returns the
    %    time vector (tout) and the output signals (yout).
    %
    % Inputs:
    %    model_name - Name of the Simulink model.
    %    T - Simulation time.
    %    x0 - Initial state vector.
    %    u - Input signal.
    %
    % Outputs:
    %    tout - Time vector.
    %    yout - Output signals.
    %
    % Example:
    %    [tout, yout] = runSimulink('Autotrans_shift', 10, [0;1000;1], [0:10]'.*ones(3,11)');
    %
    % See also: falsify
    %
    % Author: Abdelrahman Hekal
    % Written: 28-February-2023
    % Last update: ---
    % Last revision: ---

    % TODO: test and generalize function for more simulink models

    % Open the Simulink model and set the simulation mode
    open_system(model_name, 'loadonly');
    % Create Simulink.SimulationInput object (avoids permanent changes to model)
    simIn = Simulink.SimulationInput(model_name);
    
    %assign initial state
    if ~isempty(x0)
        simIn = simIn.setInitialState(x0);
    end

    % Configure input if the model has inputs
    if ~isempty(u) && ~isempty(find_system(model_name, 'BlockType', 'Inport'))
        simIn = simIn.setExternalInput(u);
    end

    % Configure simulation parameters
    assignin('base','T',T);
    simIn = simIn.setModelParameter('InitInArrayFormatMsg', 'none'); %turn off warning of initial set is array
    simIn = simIn.setModelParameter('StopTime', 'T');
    simIn = simIn.setModelParameter('SaveTime', 'on', 'TimeSaveName', 'tout');
    simIn = simIn.setModelParameter('SaveOutput', 'on', 'OutputSaveName', 'yout');
    simIn = simIn.setModelParameter('SaveFormat', 'Array');
    %save final state
    simIn = setModelParameter(simIn,"SaveFinalState","on");
    simIn = setModelParameter(simIn,"SaveOperatingPoint","on");


    % Run simulation
    simOut = sim(simIn);
    
    % Access simulation results
    tout = simOut.tout;
    yout = simOut.yout;
    %store final simout state if exists
    if isfield(simOut,'xFinal')
        xFinal=simOut.xFinal;
    else
        xFinal=[];
    end
end
