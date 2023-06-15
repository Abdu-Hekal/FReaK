function [x0,u] = getSampleXU(kfModel)
if kfModel.bestSoln.rob==inf %no previous solution, i.e. first iteration or kfModel.trainRand~=2
    [x0,u]=getRandomSampleXU(kfModel);
else
    [x0,u]=getDispSampleXU(kfModel);
end
end

function [x0,u] = getRandomSampleXU(kfModel)
%generate random initial set
x0 = randPoint(kfModel.R0);
%generate random input if kfModel has input.
cp = kfModel.cp*(kfModel.dt/kfModel.ak.dt); %get new cp array based on koopman timestep
if ~isempty(kfModel.U)
    all_steps = kfModel.T/kfModel.dt;
    u = randPoint(kfModel.U,all_steps)';
    if kfModel.pulseInput
        u = u.*kfModel.cpBool;
    else %piecewise constant input
        for k=1:length(cp)
            uk = round(cp(k)); %round is necassary for large numbers
            rep=length(u(:,k))/uk;
            assert(floor(rep)==rep,"All time steps must be a factor of cp")
            u(:,k) = repelem(u(1:uk,k),rep);
        end
    end
    u = [linspace(0,kfModel.T-kfModel.dt,all_steps)',u];
else
    u = [];
end
end

function [x0,u]=getDispSampleXU(kfModel) %TODO: implement control points
    % current values of input and initial state and valid ranges
    u = kfModel.bestSoln.u(:,2:end);
    x0 = kfModel.bestSoln.x(1,:)';
    uRange = kfModel.U;
    x0Range = kfModel.R0;
    % dimensions of u and x0
    u1 = size(u, 1);      % Number of time points
    u2 = size(u, 2);      % Number of inputs
    %displacement ratio
    dispL=0.75;
    % reshape u for size 1
    u = reshape(u,[],1);
    % Calculate normalized current sample
    curSampleOrig = [u; x0];
    uRange=repelem(uRange,u1,1); % input range for each time point
    curSample = (curSampleOrig - [uRange.inf; x0Range.inf]) ./ ([uRange.sup; x0Range.sup] - [uRange.inf; x0Range.inf]);

    % Generate random unit vector
    m = numel(u);
    n = numel(x0);
    rUnitVector = randn(m+n, 1);
    rUnitVector = rUnitVector / norm(rUnitVector);

    % Compute offsets along the random unit vector
    lam1 = ([uRange.inf; x0Range.inf] - curSample) ./ rUnitVector;
    lam2 = ([uRange.sup; x0Range.sup] - curSample) ./ rUnitVector;

    % Ensure lam1 <= lam2 for each dimension
    [lam1, lam2] = deal(min(lam1, lam2), max(lam1, lam2));
    % Handle restricted dimensions
    Ix_nonzero = find([uRange.inf; x0Range.inf] - [uRange.sup; x0Range.sup]);
    l1 = max(lam1(Ix_nonzero));
    l2 = min(lam2(Ix_nonzero));

    % Generate a random number within the range
    r = mcGenerateRandomNumber(l1, l2, dispL);

    % Create the weight vector
    weightVector = zeros(m+n, 1);
    weightVector(Ix_nonzero) = r;

    % Transform new point back to the original search hypercube
    newSample = curSample + weightVector .* rUnitVector;
    newSample = newSample .* ([uRange.sup; x0Range.sup] - [uRange.inf; x0Range.inf]) + [uRange.inf; x0Range.inf];

    % Extract the new u and newX0
    newSample(isnan(newSample)) = curSampleOrig(isnan(newSample)); %replace all nan values with original values, nan means range is zero
    newU = newSample(1:m);
    u= [kfModel.bestSoln.u(:,1),reshape(newU,u1,u2)];
    x0 = newSample(m+1:end);
end
