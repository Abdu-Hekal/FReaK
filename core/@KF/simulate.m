function [tout, yout, simTime,xFinal] = simulate(obj, x0, u, T0)
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
%    T0 (optional) - start time, default T0=0;
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
xFinal=[]; %initialize final point
if nargin<4
    T0=0;
end
sim=tic;
if isa(obj.model, 'string') || isa(obj.model,"char")
    %skip passing x0 as it is exact and set in the model. do not skip x0 if
    %it is a model operating point, it is used for partial simulations
    % TODO: check if needs to be passed
    if T0>0
        assert(isa(x0, 'Simulink.op.ModelOperatingPoint'),'if T0!=0 and model is simulink, then x0 must be a model operating point to start simulation from')
        assert(T0==x0.snapshotTime,'T0 must be equal to snapshot time of model operating point provided for starting simulation')
    end
    if ~isa(x0, 'Simulink.op.ModelOperatingPoint') && all(rad(obj.R0) == 0)
        x0=[];
    end
    [tout, yout,xFinal] = runSimulink(obj.model, obj.T, x0, u);
elseif isa(obj.model,'function_handle')
    %function handle must have 3 inputs T,x0,u
    numInputs = nargin(obj.model);
    if numInputs==2
        [tout, yout] = obj.model(obj.T-T0, x0);
    elseif numInputs==3
        [tout, yout] = obj.model(obj.T-T0, x0, u);
    else
        error(['blackbox function handle must accept 2 or 3 input ' ...
            'arguments, T, x0 and (optional) u'])
    end
    numOutputs = nargout(obj.model);
    assert(numOutputs == 2, ['blackbox function handle must return ' ...
        'two output column vectors,time points and corresponding states']);
elseif isa(obj.model,'OdeFcn')
    [tout,yout]=simulateODE(obj.model,[T0,obj.T],x0,u,obj.inputInterpolation);
elseif isa(obj.model,'ode') %object added to matlab R2023b
    F=obj.model;
    F.InitialValue = x0;
    if nargin(F.ODEFcn) > 2 %odeFcn has inputs
        F.ODEFcn = @(t, x)  F.ODEFcn(t, x, interp1(u(:,1), u(:,2:end), t, obj.inputInterpolation, 'extrap')');
    end
    sol = solve(F,T0,obj.T);
    tout=sol.Time';
    yout=sol.Solution';
else
    error('model not supported')
end
simTime=toc(sim);
%xFinal for simulink objects is a simulink model operating point used for
%partial sims, otherwise xFinal is just last yout value.
if isempty(xFinal)
    xFinal=yout(end,:)';
end
end