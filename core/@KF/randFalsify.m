function soln=randFalsify(obj)
% randFalsify - Try falsifying using random simulations for a Koopman Falsification object.
%
% Syntax:
%    soln = randFalsify(obj)
%
% Description:
%    This function attempts to falsify a set of Signal Temporal Logic (STL)
%    specifications using random simulations for a Koopman Falsification (KF)
%    object. It performs multiple simulations (up to maxSims defined), 
%    computes the robustness of the STL specifications, and returns a
%    solution struct with robustness value number of simulations, and time
%    taken
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%   soln
%       .rob - Minimum robustness value obtained during the simulations
%       .sim     - Number of simulations completed
%       .runtime - time taken
%
% Example:
%    soln = randFalsify(obj);
%
% See also: randSimulation, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

tic;
soln.falsified=false;
soln.rob=inf;
tsim = (0:obj.dt:obj.T)'; %define time points for interpolating simulation
for sim=1:obj.maxSims
    [t, x, ~, u] = randSimulation(obj);
    if ~isempty(u)
        interpU = interp1(u(:,1),u(:,2:end),t,obj.inputInterpolation); %interpolate input at same time points as trajectory
    else
        interpU=u;
    end

    for jj=1:length(obj.spec)
        [~,~,robustness] = bReachRob(obj.spec(jj),tsim,x,interpU(:,2:end)');
        if robustness < soln.rob
            soln.rob=robustness;
        end
        if soln.rob<0
            soln.falsified=true;
            break;
        end
    end
end
soln.sims=sim; soln.runtime=toc;
LogicalStr = {'No', 'Yes'};
vprintf(obj.verb,1,'<-------------------------------------------------------> \n')
vprintf(obj.verb,1,"Falsified: %s \n",LogicalStr{soln.falsified+1})
vprintf(obj.verb,1,'Minimum robustness: %.2f after %d simulations \n',soln.rob,soln.sims);
vprintf(obj.verb,1,'Time taken: %.2f seconds\n',soln.runtime);
end