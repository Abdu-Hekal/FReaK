function [tout, yout, te, ye, ie] = simulateODE(obj, tspan, x0, u, inputInterpolation)
    % simulateODE - Simulate ODEs using the provided OdeFcn object
    %
    % Syntax:
    %    [tout, yout] = simulateODE(obj, tspan, x0, u)
    %    (additional outputs for events):
    %    [tout, yout, te, ye, ie] = simulateODE(obj, tspan, x0, u)
    %
    % Inputs:
    %    obj - OdeFcn object
    %    tspan - Time span for simulation [t_start, t_end]
    %    x0 - Initial conditions for the ODEs
    %    u - Matrix of inputs where the first column is time and subsequent columns are inputs
    %
    % Outputs:
    %    tout - Output time vector of the simulation
    %    yout - Output state vector of the simulation
    %    (The following outputs are in case there is a detected event, 
    %       see https://uk.mathworks.com/help/matlab/math/ode-event-location.html)
    %    te - a column vector of the times at which events occurred.
    %    ye - the solution value at each of the event times in te.
    %    ie - contains indices into the vector returned by the event 
    %         function. The values indicate which event the solver detected.
    %
    % Example:
    % odeFun = @(t, x, u) -x + u; % Example ODE function with input
    % solver = 'ode45';
    % options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    % odeSim = OdeFcn(odeFun, solver, options);
    % tspan = [0, 10];
    % x0 = 1;
    % t_values=[0;1];
    % input_values=[0;1];
    % u = [t_values, input_values]; % Input matrix
    % [tout, yout] = simulateODE(odeSim, tspan, x0, u);
    % 
    % See also: OdeFcn, ode45, ode23, odeset
    % Extract handle from the OdeFcn object
    odeFcn = obj.handle;

    % Define the ODE function with input
    if nargin(odeFcn) > 2 %odeFcn has inputs
        assert(size(u,1)>=2,'Input must have at least two sample points')
        assert(size(u,2)>=2,'Input must have at least two columns, where first column is time points')
        if nargin < 5 %no interpolation scheme defined
            inputInterpolation='pchip';
        end
        odeFcn = @(t, x) odeFcn(t, x, interp1(u(:,1), u(:,2:end), t, inputInterpolation, 'extrap')');
    end

    % Perform ODE simulation
    if nargout>2
        assert(~isempty(obj.options.Events),'ODE must have events to return additional event outputs')
        [tout, yout, te, ye, ie] = feval(obj.solver, odeFcn, tspan, x0, obj.options);
    else
        [tout,yout] = feval(obj.solver, odeFcn, tspan, x0, obj.options);
    end
end
