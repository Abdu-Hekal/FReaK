function [tout, yout, simTime] = simulate(obj, x0, u)
% simulate - Simulate the model associated with a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, simTime] = simulate(obj, x0, u)
%
% Description:
%    This function simulates the model associated with a Koopman
%    Falsification (KF) object. The model can be either a Simulink model
%    or a custom function handle. Note that a custom function can be
%    used for simulation by passing the function handle. Ensure that the
%    outputs are consistent with the simulation method used for the
%    specific model.
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
% See also: falsify, sampleSimulation
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------
tic;
if isa(obj.model, 'string') || isa(obj.model,"char")
    %skip passing x0 as it is exact and set in the model. TODO: check if
    %needs to be passed
    if all(rad(obj.R0) == 0)
        x0=[];
    end
    [tout, yout] = runSimulink(obj.model, obj.T, x0, u);
elseif isa(obj.model,'function_handle')
    %function handle must have 3 inputs T,x0,u
    numInputs = nargin(obj.model);
    if numInputs==2
        [tout, yout] = obj.model(obj.T, x0);
    elseif numInputs==3
        [tout, yout] = obj.model(obj.T, x0, u);
    else
        error(['blackbox function handle must accept 2 or 3 input ' ...
            'arguments, T, x0 and (optional) u'])
    end
    numOutputs = nargout(obj.model);
    assert(numOutputs == 2, ['blackbox function handle must return ' ...
        'two output column vectors,time points and corresponding states']);
elseif isa(obj.model,'OdeFcn')
    [tout,yout]=simulateODE(obj.model,[0,obj.T],x0,u,obj.inputInterpolation);
elseif isa(obj.model,'ode') %object added to matlab R2023b
    F=obj.model;
    F.InitialValue = x0;
    if nargin(F.ODEFcn) > 2 %odeFcn has inputs
        F.ODEFcn = @(t, x)  F.ODEFcn(t, x, interp1(u(:,1), u(:,2:end), t, obj.inputInterpolation, 'extrap')');
    end
    sol = solve(F,0,obj.T);
    tout=sol.Time';
    yout=sol.Solution';
else
    error('model not supported')
end
simTime=toc;
end