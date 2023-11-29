function [tout, yout] = runSimulink(model_name, T, x0, u)
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

%TODO: test and generalize function for more simulink models
%------------- BEGIN CODE --------------

open_system(model_name,'loadonly');
sim_mode = get_param(model_name,'SimulationMode');
set_param(model_name,'SimulationMode','normal'); %need this to extract init state info, otherwise error is thrown if sim mode is not 'normal'

xInitial = Simulink.BlockDiagram.getInitialState(model_name);
for i=1:size(xInitial.signals,2)
    xInitial.signals(i).values = x0(i);
end
paramStruct.LoadInitialState = 'on';
paramStruct.InitialState = 'xInitial';
assignin('base','T',T);
assignin('base','xInitial',xInitial);

if ~isempty(find_system(model_name,'BlockType','Inport')) %check if simulink model has inputs
    paramStruct.LoadExternalInput= 'on';
    paramStruct.ExternalInput= 'u';
    assignin('base','u',u);
end

paramStruct.StopTime = 'T';
paramStruct.SaveTime='on';
paramStruct.TimeSaveName='tout';
paramStruct.SaveOutput='on';
paramStruct.OutputSaveName='yout';
paramStruct.SaveFormat='Array';
paramStruct.SimulationMode=sim_mode;

result = sim(model_name, paramStruct);
tout = result.tout;
yout = result.yout;

end
