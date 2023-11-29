function [minRob,ii]=randFalsify(obj)
% randFalsify - Try falsifying using random simulations for a Koopman Falsification object.
%
% Syntax:
%    [minRob, ii] = randFalsify(obj)
%
% Description:
%    This function attempts to falsify a set of Signal Temporal Logic (STL)
%    specifications using random simulations for a Koopman Falsification (KF)
%    object. It performs multiple simulations (up to maxSims defined), 
%    computes the robustness of the STL specifications, and returns the
%    minimum robustness value along with the corresponding number of simulations.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%    minRob - Minimum robustness value obtained during the simulations
%    ii     - Simulation index corresponding to the minimum robustness
%
% Example:
%    [minRob, ii] = randFalsify(obj);
%
% See also: randSimulation, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

minRob=inf;
tsim = (0:obj.dt:obj.T)'; %define time points for interpolating simulation
for ii=1:obj.maxSims
    [t, x, ~, u] = randSimulation(obj);
    if ~isempty(u)
        interpU = interp1(u(:,1),u(:,2:end),t,obj.inputInterpolation); %interpolate input at same time points as trajectory
    else
        interpU=u;
    end

    for jj=1:length(obj.spec)
        [~,~,robustness] = bReachRob(obj.spec(jj),tsim,x,interpU(:,2:end)');
        if robustness < minRob
            minRob=robustness;
        end
        if minRob<0
            fprintf('Falsified after %d iterations \n',ii);
            fprintf('Minimum robustness: %f \n',minRob);
            return;
        end
    end
end
fprintf('Failed to falsify in %d iterations \n',ii);
fprintf('Minimum robustness: %f \n',minRob);
end