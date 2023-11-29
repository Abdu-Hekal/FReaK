function [tout, yout, x0, u, obj] = randSimulation(obj)
% randSimulation - Perform a random simulation for a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, x0, u, obj] = randSimulation(obj)
%
% Description:
%    This function generates a random initial state and input, then performs
%    a simulation for a Koopman Falsification (KF) object. The simulated
%    time vector, output vector, updated initial state, input vector, and
%    the updated KF object are returned.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    x0   - Updated initial state vector after simulation
%    u    - Input vector used for simulation
%    obj  - Updated Koopman Falsification (KF) object with simulation
%           statistics.
%
% Example:
%    [tout, yout, x0, u, obj] = randSimulation(obj);
%
% See also: simulate, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

[x0,u] = getRandomSampleXU(obj);
tsim = (0:obj.dt:obj.T)'; %define time points for interpolating input
if ~isempty(u)
    usim = interp1(u(:,1),u(:,2:end),tsim,obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
    usim =  max(obj.U.inf',min(obj.U.sup',usim)); %ensure that extrapolation is within input bounds
    usim = [tsim,usim];
else
    usim=u; %no input for the model
end
[tout, yout, obj] = simulate(obj, x0, usim);
end
