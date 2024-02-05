function [tout, yout, u, simTime] = sampleSimulation(obj,varargin)
% sampleSimulation - Perform a random simulation for a Koopman Falsification object.
%
% Syntax:
%    [tout, yout, u, simTime] = sampleSimulation(obj)
%    [tout, yout, u, simTime] = sampleSimulation(obj,allData)
%    [tout, yout, u, simTime] = sampleSimulation(obj,allData,perturb)
%
% Description:
%    This function generates a sample initial state and input, then performs
%    a simulation for a Koopman Falsification (KF) object. The simulated
%    time vector, output vector, updated initial state, input vector, and
%    the time taken for simulation are returned.
%   The sample is either random, a perturbation of the current best
%   sample or new sample using mopso (from staliro soar) depending on
%   KF settings and passed arguments
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

if isempty(obj.inputsInterval)
    obj=initialize(obj);
    vprintf(obj.verb,1,'Initiallizing object to setup control points and input intervals\n')
end

perturb=0;
simd=true; %simulation needed
rnd=true; %find a random point if needed
if obj.resetStrat>=1
    try
        if nargin>1
            allData=varargin{1};
            assert(isstruct(allData), 'all Data must be a struct')
        end
        if nargin>2
            assert(isnumeric(perturb), 'perturb must be a numeric value');
            assert(perturb >= 0 && perturb <= 1, 'perturb must be between 0 and 1');
            perturb=varargin{2};
        end

        [bestRob,bestIdx] = min(allData.Rob);
        if ~(isempty(allData.X)||bestRob==inf) %if no input samples then random sampling is used
            allX0=cellfun(@(x) x(1, :), allData.X, 'UniformOutput', false);
            allInputsSamples=cell2mat(allX0');
            if ~isempty(obj.U)
                inputs = cellfun(@(u) reshape(u(1:end-1,2:end), [], 1), allData.XU, 'UniformOutput', false);
                inputs = cellfun(@(u) u(find(obj.cpBool)),inputs,'UniformOutput', false);
                inputs = cell2mat(inputs)';
                allInputsSamples=[allInputsSamples,inputs];
            end
            if obj.resetStrat==1
                if perturb==0 %return exact bestSoln as new training sample
                    tout=allData.t{bestIdx};
                    yout=allData.X{bestIdx};
                    u=allData.XU{bestIdx};
                    simTime=0;
                    simd=false; %no simulation needed
                else
                    bestSample=allInputsSamples(bestIdx,:)';
                    sample=getDispSampleXU(obj,bestSample,perturb);
                end
                rnd=false;
            end
            if obj.resetStrat==2
                sample=getMopsoSampleXU(obj,allInputsSamples,allData.Rob);
                rnd=false;
            end
        end
    catch ME
        vprintf(obj.verb,2,'Error in reset strategy %d, due to "%s", using random reset \n',obj.resetStrat,ME.message)
    end
end

if rnd %default random simulations
    sample = randPoint(obj.inputsInterval)';
end
if simd %simulate sample if needed
    [x0,u]=getInputs(obj,sample);
    [tout, yout, simTime] = simulate(obj, x0, u);
end
end

function newSample=getDispSampleXU(obj,bestSample,perturb)

lowerBound = obj.inputsInterval.inf;
upperBound = obj.inputsInterval.sup;

maxPerturb = perturb * (upperBound-lowerBound);

lowerBound = max(bestSample-maxPerturb,lowerBound); %maximum of perturbation and bounds on inputs
upperBound = min(bestSample+maxPerturb,upperBound); %minimum of perturbation and bounds on inputs

newSample = (upperBound - lowerBound) .* rand(size(bestSample)) + lowerBound;
end

function sample=getMopsoSampleXU(obj,allInputsSamples,allRob)

%find non-exact dimensions as covariance does not perform well for exact dims
nonExactDims=find(rad(obj.inputsInterval));
allInputsSamples=allInputsSamples(:,nonExactDims);


% Initialize MOPSO Parameters
MOparams.Np = 200;        % Population size
MOparams.Nr = 200;        % Repository size
MOparams.maxgen = 500;    % Maximum number of generations
MOparams.W = 0.4;         % Inertia weight
MOparams.C1 = 2;          % Individual confidence factor
MOparams.C2 = 2;          % Swarm confidence factor
MOparams.ngrid = 20;      % Number of grids in each dimension
MOparams.maxvel = 5;      % Maxmium vel in percentage
MOparams.u_mut = 0.5;     % Uniform mutation percentage
MultiObj.nVar = size(nonExactDims,1);  % Set problem dimension
MultiObj.var_min = obj.inputsInterval.inf(nonExactDims)';
MultiObj.var_max = obj.inputsInterval.sup(nonExactDims)';


%Fit Gaussian Process Meta Mod  el
GPmod = OK_Rmodel_kd_nugget(allInputsSamples, allRob, 0, 2);

% optimize EI and CD with MOPSO
MultiObj.fun = @(x)[-EIcalc_kd(x,allInputsSamples,GPmod,allRob), -CrowdingDist_kd(x,allInputsSamples)];

pf = MOPSO(MOparams,MultiObj);
[minNegEI, index] = min(pf.pos_fit(:,1));

crowded_EI_flag=1;
alpha_lvl_set = 0.05; %the top percent of EI values to be inluded in the crowding distance phase
%use crowded EI
if crowded_EI_flag == 1
    best_crowd = inf;
    for k = 1:size(pf.pos,1)
        if pf.pos_fit(k,1) <= (minNegEI*(1-alpha_lvl_set))
            if pf.pos_fit(k,2) < best_crowd
                best_crowd = pf.pos_fit(k,2);
                newSample = pf.pos(k,:);
            end
        end
    end
    %use standard EI
else
    newSample = pf.pos(index,:);
end
sample=center(obj.inputsInterval)';
sample(nonExactDims)=newSample;
end

function [x0,u]=getInputs(obj,sample)
x0 = sample(1:dim(obj.R0))';
if ~isempty(obj.U)
    u=[];
    allSteps = obj.T/obj.ak.dt;
    allTimePoints=linspace(0,obj.T,allSteps+1)';

    newU = sample(dim(obj.R0)+1:end);
    assert(numel(newU)==sum(obj.cp),'Number of inputs does not match with generated samples, check for errors')
    
    for i=1:length(obj.cp)
        cp=obj.cp(:,i);
        if cp==1 %only one control input, apply to beginning an end
            timePoints=[0;obj.T];
            cpVal=[newU(1);newU(1)];
        else
            timePoints=linspace(0,obj.T,cp+1)';
            timePoints = timePoints(1:end-1); %remove last time point, inputs are up to step k-1
            cpVal=newU(1:cp);
        end
        u(:,i) = interp1(timePoints, cpVal, allTimePoints,obj.inputInterpolation,"extrap");
        newU = newU(cp+1:end); %skip to next input dimension
    end
    u= [(0:obj.ak.dt:obj.T)',u]; %append time
else
    u=[];
end
end

