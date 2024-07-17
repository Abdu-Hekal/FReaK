function varargout = randFalsify(obj)
% randFalsify - Try falsifying using random simulations for a Koopman Falsification object.
%
% Syntax:
%    soln = randFalsify(obj)
%    [soln,sims] = randFalsify(obj)
%
% Description:
%    This function attempts to falsify a set of Signal Temporal Logic (STL)
%    specifications using random simulations for a Koopman Falsification (KF)
%    object. It performs multiple simulations (up to maxSims defined), 
%    computes the robustness of the STL specifications, and returns a
%    solution struct with robustness value number of simulations, and time
%    taken. It also optionally returns a struct sims with data of all
%    simulations completed and corresponding robustness
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%   soln
%       .best - struct with best soln found
%           .rob - Minimum robustness value 
%           .t - array of time points
%           .x - trajectory
%           .u - array of inputs
%
%       .sim     - Number of simulations completed
%       .runtime - time taken
%   sims
%       .t - set of time points
%       .X - set of trajectories
%       .XU - set of inputs
%       .ROB - set of robustness values
%
% Example:
%    soln = randFalsify(obj);
%
% See also: sampleSimulation, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

tic;
soln.falsified=false;
soln.best.rob=inf;
soln.sims=0;
if nargout > 1
    sims.t={}; sims.X={}; sims.XU={}; sims.ROB={};
end
for ii=1:obj.maxSims
    [t, x, u] = sampleSimulation(obj);
    soln.sims = soln.sims+1; 
    [soln,falsified,robustness]=checkFalsification(soln,x,u,t,obj.spec,obj.inputInterpolation,'rand Simulation',obj.verb);
    if robustness < soln.best.rob
        soln.best.rob=robustness; soln.best.t=t; soln.best.x=x; soln.best.u=u;
    end
    if nargout > 1
        sims.t{end+1}=t; sims.X{end+1}=x; sims.XU{end+1}=u; sims.ROB{end+1}=robustness;
    end
    if falsified
        soln.falsified=true;
        break;
    end

end
if nargout > 0
    varargout{1} = soln;
end
if nargout > 1
    varargout{2} = sims;
end

soln.runtime=toc;
LogicalStr = {'No', 'Yes'};
vprintf(obj.verb,1,'<-------------------------------------------------------> \n')
vprintf(obj.verb,1,"Falsified: %s \n",LogicalStr{soln.falsified+1})
vprintf(obj.verb,1,'Minimum robustness: %.2f after %d simulations \n',soln.best.rob,soln.sims);
vprintf(obj.verb,1,'Time taken: %.2f seconds\n',soln.runtime);
end