function [soln,falsified,robustness,Bdata,newGlobalBest,bestSpec]=checkFalsification(soln,x,u,t,specs,tcp,inputInterpolation,U,method,verb)
% CHECKFALSIFICATION Checks if a given trajectory falsifies a set of specifications.
%
%   [soln, falsified, robustness, Bdata, newGlobalBest] = checkFalsification(soln, x, u, t, specs, inputInterpolation, method, verb)
%
% Inputs:
%   soln - Struct containing information about the current falsification process.
%   x - State trajectory.
%   u - Input trajectory.
%   t - Time points corresponding to the trajectories.
%   specs - CORA specification object.
%   inputInterpolation - Interpolation method for input trajectories.
%   method - String indicating the method used for falsification.
%   verb - Verbosity level for printing messages.
%
% Outputs:
%   soln - Updated struct with information about the current falsification process.
%   falsified - Logical indicating if the trajectory falsifies any specification.
%   robustness - Robustness value associated with the falsification.
%   Bdata - Data associated with the falsification, e.g., predicates and times.
%   newGlobalBest - Logical indicating if a new best solution  has been found.
%
% Description:
%   This function checks if a given trajectory falsifies a set of specifications.
%   The specifications can be of different types, including unsafe sets, safe sets,
%   and logical formulas. The function evaluates the trajectory against each
%   specification and updates the falsification status, robustness value,
%   and relevant data. Additionally, it checks if the current solution is better
%   than the best solution found so far and updates the soln struct accordingly.
%
% See also:
%   bReachRob, interp1, vprintf
%
% Author: Abdelrahman Hekal
% Written: 28-February-2023
% Last update: ---

falsified=false; bestRob=inf; bestSpec=[];
Bdata=NaN; newGlobalBest=false;
for ii=1:numel(specs)
    spec=specs(ii);
    % different types of specifications
    %TODO, compute robustness for sets
    if strcmp(spec.type,'unsafeSet')
        falsified = any(spec.set.contains(x'));
        robustness=inf;
    elseif strcmp(spec.type,'safeSet')
        falsified = ~all(spec.set.contains(x'));
        robustness=inf;
    elseif strcmp(spec.type,'logic')
        if ~isempty(u)
            interpU = interp1(u(:,1),u(:,2:end),t,inputInterpolation); %interpolate input at same time points as trajectory
            usim = interp1(u(:,1),u(:,2:end),tcp,inputInterpolation,"extrap"); %interpolate and extrapolate input points
            usim = max(U.inf',min(U.sup',usim)); %ensure that extrapolation is within input bounds
        else
            interpU=u;
            usim=u;
        end
        [Bdata,~,robustness] = bReachRob(spec,t,x,interpU');
        vprintf(verb,3,"robustness value: %.3f after %d simulations with method: %s \n",robustness,soln.sims,method)
        if robustness <= 0 %if less than or equal zero
            falsified=true;
        end
    end
    if robustness<=bestRob
        bestRob=robustness;
        bestSpec=spec;
    end
    if robustness < soln.best.rob
        vprintf(verb,2,"new best robustness!: %.3f after %d simulations due to: %s \n",robustness,soln.sims,method)
        soln.best.rob=robustness;
        soln.best.x=x; soln.best.u=usim; soln.best.t=tcp; soln.best.spec=spec;
        newGlobalBest=true; %found a new best soln
    end
    if falsified
        break;
    end
end
end