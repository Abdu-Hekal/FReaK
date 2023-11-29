function [tout, yout, obj] = simulate(obj, x0, u)
% simulate - Simulate the model associated with a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, obj] = simulate(obj, x0, u)
%
% Description:
%    This function simulates the model associated with a Koopman
%    Falsification (KF) object. The model can be either a Simulink model
%    or a custom function handle. Note that a custom simulate function can be
%    used for a subclass of this class. Ensure that the outputs are consistent
%    with the simulation method used for the specific model.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%    x0  - Initial state vector for simulation
%    u   - Input vector for simulation
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    obj  - Updated Koopman Falsification (KF) object with simulation
%           statistics.
%
% Example:
%    [tout, yout, obj] = simulate(obj, x0, u);
%
% See also: falsify, randSimulation
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

tic
if isa(obj.model, 'string') || isa(obj.model,"char")
    [tout, yout] = runSimulink(obj.model, obj.T, x0, u);
elseif isa(obj.model,'function_handle')
    %function handle must have 3 inputs T,x0,u
    [tout, yout] = obj.model(obj.T, x0, u);
else
    error('model not supported')
end
obj.soln.sims = obj.soln.sims+1;
obj.soln.simTime = obj.soln.simTime+toc;
end