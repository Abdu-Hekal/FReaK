function [x0,u] = getSampleXU(kfModel)
if kfModel.bestSoln.rob==inf %no previous solution, i.e. first iteration or kfModel.trainRand~=2
    [x0,u]=getRandomSampleXU(kfModel);
else
    [x0,u]=getDispSampleXU(kfModel);
end
end

function [x0,u]=getDispSampleXU(kfModel) %TODO: implement control points

u = kfModel.bestSoln.u(:,2:end);
x0 = kfModel.bestSoln.x(1,:)';
uRange = kfModel.U;
x0Range = kfModel.R0;
u1 = size(u, 1);      % Number of time points
u2 = size(u, 2);      % Number of inputs

perturb = 0.1; %max perturbation percentage
lowerBound = [repmat(uRange.inf,size(u,1),1); x0Range.inf];
upperBound = [repmat(uRange.sup,size(u,1),1); x0Range.sup];

u = reshape(u,[],1);
curSample = [u; x0];

maxPerturb = perturb * (upperBound-lowerBound);
lowerBound = max(curSample-maxPerturb,lowerBound); %maximum of perturbation and bounds on inputs
upperBound = min(curSample+maxPerturb,upperBound); %minimum of perturbation and bounds on inputs

newSample = (upperBound - lowerBound) .* rand(size(curSample)) + lowerBound;
newU = newSample(1:u1*u2);
u= [kfModel.bestSoln.u(:,1),reshape(newU,u1,u2)]; %append time from previously found best solution
x0 = newSample(u1*u2+1:end);

end
