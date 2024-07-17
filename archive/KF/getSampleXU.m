function [x0,u] = getSampleXU(obj,bestSoln)
% getSampleXU - Generate a sample of initial state (x0) and input signal (u)
%
% Syntax:
%    [x0, u] = getSampleXU(obj,bestSoln)
%
% Description:
%    This function generates a sample of initial state (x0) and input signal (u)
%    based on the current state of the Koopman Falsification (KF) model. If
%    the KF model has a previous best solution, it uses a perturbed version
%    of that solution. Otherwise, it generates a random sample.
%
% Inputs:
%    obj - KF object containing the Koopman model and various parameters
%              needed for the falsification process.
%   bestSoln - bestSoln found so far
%
% Outputs:
%    x0 - Sampled initial state.
%    u - Sampled input signal.
%
% Example:
%    [x0, u] = getSampleXU(obj);
%
% See also: getRandomSampleXU
%
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: [Date]
% Last revision: [Date]
%
% ------------- BEGIN CODE --------------

if obj.trainStrat~=2 || bestSoln.rob==inf %no previous solution, i.e. obj.trainStrat~=2 or first iteration
    [x0,u]=getRandomSampleXU(obj);
else
    [x0,u]=getDispSampleXU(obj,bestSoln);
end
end

function [x0,u]=getDispSampleXU(obj,bestSoln) %TODO: implement control points

u = bestSoln.u(:,2:end);
x0 = bestSoln.x(1,:)';
uRange = obj.U;
x0Range = obj.R0;
u1 = size(u, 1);      % Number of time points
u2 = size(u, 2);      % Number of inputs

perturb = 0.01; %max perturbation percentage
lowerBound = [repmat(uRange.inf,size(u,1),1); x0Range.inf];
upperBound = [repmat(uRange.sup,size(u,1),1); x0Range.sup];

u = reshape(u,[],1);
curSample = [u; x0];

maxPerturb = perturb * (upperBound-lowerBound);

lowerBound = max(curSample-maxPerturb,lowerBound); %maximum of perturbation and bounds on inputs
upperBound = min(curSample+maxPerturb,upperBound); %minimum of perturbation and bounds on inputs

newSample = (upperBound - lowerBound) .* rand(size(curSample)) + lowerBound;
newU = newSample(1:u1*u2);
u= [bestSoln.u(:,1),reshape(newU,u1,u2)]; %append time from previously found best solution
x0 = newSample(u1*u2+1:end);

end
