function [tout, yout, u, simTime,perturb] = sampleSimulation(obj,varargin)
% sampleSimulation - Perform a random simulation for a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, u, simTime] = sampleSimulation(obj)
%
% Description:
%    This function generates a sample initial state and input, then performs
%    a simulation for a Koopman Falsification (KF) object. The simulated
%    time vector, output vector, updated initial state, input vector, and
%    the time taken for simulation are returned.
%   The sample is either random or a perturbation of the current best
%   sample depending on KF settings and passed arguments
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%    bestSoln (optional) - best soln found so far. used for generating
%    perturbed sample
%    perturb (optional) - perturbation factor (default=0)
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    u    - Input vector used for simulation
%    simTime  - time taken for simulation
%
% Example:
%    [tout, yout, u, simTime] = sampleSimulation(obj);
%
% See also: simulate, falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

bestSoln=[]; perturb=0;
if nargin>1
    bestSoln=varargin{1};
    assert(isstruct(bestSoln), 'provided argument must be a struct storing best solution');
    % Check if the struct has fields x, u, and t
    assert(isfield(bestSoln, 'x'), 'bestSoln must have a field named x');
    assert(isfield(bestSoln, 'u'), 'bestSoln must have a field named u');
    assert(isfield(bestSoln, 't'), 'bestSoln must have a field named t');
end
if nargin>2
    assert(isnumeric(perturb), 'perturb must be a numeric value');
    assert(perturb >= 0 && perturb <= 1, 'perturb must be between 0 and 1');
    perturb=varargin{2};
end
if isempty(bestSoln) || obj.trainStrat==0 || obj.trainStrat==3 || bestSoln.rob==inf
    [x0,u] = getRandomSampleXU(obj);
    [tout, yout, simTime] = simulate(obj, x0, u);
else
    if perturb==0 %return exact bestSoln as new training sample
        tout=bestSoln.t;
        yout=bestSoln.x;
        u=bestSoln.u;
        simTime=0;
    else
        [x0,u]=getDispSampleXU(obj,bestSoln,perturb);
        [tout, yout, simTime] = simulate(obj, x0, u);
    end
    perturb=min(1,perturb+obj.sampPerturb);
end
end

function [x0,u]=getDispSampleXU(obj,bestSoln,perturb)

u = bestSoln.u(:,2:end);
x0 = bestSoln.x(1,:)';
uRange = obj.U;
x0Range = obj.R0;
u1 = size(u, 1);      % Number of time points
u2 = size(u, 2);      % Number of inputs

lowerBound = [];
upperBound=[];
for i=1:size(uRange,1)
    lowerBound=[lowerBound;repmat(uRange.inf(i),size(u,1),1)];
    upperBound = [upperBound;repmat(uRange.sup(i),size(u,1),1)];
end
lowerBound=[lowerBound;x0Range.inf];
upperBound=[upperBound;x0Range.sup];

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

