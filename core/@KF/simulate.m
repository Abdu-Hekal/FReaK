function [tout, yout, simTime] = simulate(obj, x0, u)
% simulate - Simulate the model associated with a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, simTime] = simulate(obj, x0, u)
%
% Description:
%    This function simulates the model associated with a Koopman
%    Falsification (KF) object. The model can be either a Simulink model
%    or a custom function handle. Note that a custom simulate function can be
%    used for a subclass of this class. Ensure that the outputs are consistent
%    with the simulation method used for the specific model.
%
% Inputs:
%    x0  - Initial state vector for simulation
%    u   - Input vector for simulation
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    simTime  - time taken for simulation
%
% Example:
%    [tout, yout, simTime] = simulate(obj, x0, u);
%
% See also: falsify, randSimulation
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

%interpolate input in accordance with interpolation strategy defined
tsim = (0:obj.dt:obj.T)'; %time points for interpolating input
if ~isempty(u)
    usim = interp1(u(:,1),u(:,2:end),tsim,obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
    usim =  max(obj.U.inf',min(obj.U.sup',usim)); %ensure that extrapolation is within input bounds
    usim = [tsim,usim];
else
    usim=u; %no input for the model
end

tic
if isa(obj.model, 'string') || isa(obj.model,"char")
     %skip passing x0 as it is exact and set in the model. TODO: check if
     %needs to be passed
    if all(rad(obj.R0) == 0)
        x0=[];
    end
    [tout, yout] = runSimulink(obj.model, obj.T, x0, usim);
elseif isa(obj.model,'function_handle')
    %function handle must have 3 inputs T,x0,u
    [tout, yout] = obj.model(obj.T, x0, usim);
else
    error('model not supported')
end
simTime=toc;
end