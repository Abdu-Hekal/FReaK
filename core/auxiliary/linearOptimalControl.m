function critU = linearOptimalControl(obj, koopModel, xTarget)
% LINEAROPTIMALCONTROL Computes linear optimal control inputs for tracking a target trajectory.
%
%   critU = linearOptimalControl(obj, koopModel, xTarget)
%
% Inputs:
%   obj - KF object containing system parameters and specifications.
%   koopModel - Koopman model struct with matrices A and B.
%   xTarget - Target trajectory to be tracked.
%
% Outputs:
%   critU - Optimal control inputs for tracking the target trajectory.
%
% Description:
%   This function calculates the optimal control inputs using a linear model
%   (Koopman model) to track a target trajectory. The optimization problem
%   minimizes the tracking error by adjusting the control inputs. The target
%   trajectory is filled in by linear extrapolation, and the optimization
%   considers constraints on inputs and the linear system dynamics.
%
%   The optimization problem is formulated using YALMIP and solved using the
%   specified solver in the KF object. The resulting optimal control inputs
%   are then extracted and returned. Additionally, the last time point is
%   extrapolated to cover the entire simulation duration, and the inputs are
%   interpolated for smoother trajectories.
%
% See also:
%   KF, sdpvar, interp1, optimize, value
%
% Author: Abdelrahman Hekal
% Written: 28-February-2023
% Last update: ---

% fill in xTarget by linear extrapolation, when solver step is not for
% every koopman step, xTarget is only for solver steps and rest is NaN
xTarget = arrayfun(@(row) interp1(find(~isnan(xTarget(row, :))), xTarget(row, ~isnan(xTarget(row, :))), 1:size(xTarget, 2), 'linear', 'extrap'), 1:size(xTarget, 1), 'UniformOutput', false)';
xTarget = cell2mat(xTarget);

% Define the state and control variables
n = size(xTarget, 1);     % Number of state variables
m = size(xTarget, 2);     % Number of discrete time points
nu = dim(obj.U); %number of inputs

nObs = size(koopModel.A, 1); % Number of state variables + koopman observables
x = sdpvar(nObs, m);
u = sdpvar(nu, m-1);

% Tracking error objective as squared sum of distance, TODO: add state cost
% matrix 'Q' to control importance of tracking error in each dimenstion,
% focus more on dimensions in stl.
tracking_error = x(1:n, 1:m) - xTarget(:, 1:m);
cost = sum(sum(tracking_error.^2));

% Tracking error objective using Frobenius norm
% tracking_error = x(1:n, 1:m) - xTarget(:, 1:m);
% cost = norm(tracking_error, 'fro')^2;

% Vectorized constraints for inputs and system dynamics
constraints = [u >= repmat(obj.U.inf,1,m-1), u <= repmat(obj.U.sup,1,m-1)];
constraints = [constraints, x(:, 2:m) == koopModel.A * x(:, 1:(m - 1)) + koopModel.B * u(:, 1:(m - 1))];

% Define the optimization problem
optimize(constraints, cost, obj.solver.opts);

% Extract and return the optimal control inputs
critU = value(u);
% extrapolate last time point (dummy input) and add time steps
all_steps = obj.T/obj.ak.dt;
tp_=linspace(0,obj.T-obj.ak.dt,all_steps); %time points without last time step
tp = linspace(0,obj.T,all_steps+1);
critU = interp1(tp_',critU',tp',obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
critU =  max(obj.U.inf',min(obj.U.sup',critU)); %ensure that extrapolation is within input bounds
critU = [tp',critU];

end
