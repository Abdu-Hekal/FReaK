function [tout, yout, x0, u, simTime] = randSimulation(obj)
% randSimulation - Perform a random simulation for a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, x0, u, simTime] = randSimulation(obj)
%
% Description:
%    This function generates a random initial state and input, then performs
%    a simulation for a Koopman Falsification (KF) object. The simulated
%    time vector, output vector, updated initial state, input vector, and
%    the time taken for simulation are returned.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    x0   - Updated initial state vector after simulation
%    u    - Input vector used for simulation
%    simTime  - time taken for simulation
%
% Example:
%    [tout, yout, x0, u, simTime] = randSimulation(obj);
%
% See also: simulate, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

[x0,u] = getRandomSampleXU(obj);
[tout, yout, simTime] = simulate(obj, x0, u);
end
