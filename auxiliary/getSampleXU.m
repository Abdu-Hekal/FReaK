function [x0,u] = getSampleXU(model)
if model.bestSoln.rob==inf %no previous solution, i.e. first iteration
    [x0,u]=getRandomSampleXU(model);
else
    [x0,u]=getDispSampleXU(model);
end
end

function [x0,u] = getRandomSampleXU(model)
%generate random initial set
x0 = randPoint(model.R0);
%generate random input if model has input.
if ~isempty(model.U)
    all_steps = model.T/model.dt;
    u = randPoint(model.U,all_steps)';
    if model.pulseInput
        u = u.*model.cpBool;
    else %piecewise constant input
        for k=1:length(model.cp)
            uk = model.cp(k);
            u(:,k) = repelem(u(1:uk,k),length(u(:,k))/uk);
        end
    end
    u = [linspace(0,model.T-model.dt,all_steps)',u];
else
    u = [];
end
end

function [x0,u]=getDispSampleXU(model)
    % current values of input and initial state and valid ranges
    u = model.bestSoln.u(:,2:end);
    x0 = model.bestSoln.x(1,:)';
    uRange = model.U;
    x0Range = model.R0;
    % dimensions of u and x0
    u1 = size(u, 1);      % Number of time points
    u2 = size(u, 2);      % Number of inputs
    %displacement ratio
    dispL=10.75;
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
    u= [model.bestSoln.u(:,1),reshape(newU,u1,u2)];
    x0 = newSample(m+1:end);
end
