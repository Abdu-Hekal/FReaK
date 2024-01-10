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

allU=[]; allRob=[];
if nargin>1
    allU=varargin{1};
end
if nargin>2
    allRob=varargin{2};
end
if isempty(allU) || obj.trainStrat==0 || obj.trainStrat==3
    [x0,u] = getRandomSampleXU(obj);
    [tout, yout, simTime] = simulate(obj, x0, u);
else
    try
        [x0,u]=getDispSampleXU(obj,allU,allRob);
    catch
        vprintf(obj.verb,1,'Error in mopso optimization for reset, using random reset \n')
        [x0,u] = getRandomSampleXU(obj);
    end
    [tout, yout, simTime] = simulate(obj, x0, u);
end
end

function [x0,u]=getDispSampleXU(obj,allU,allRob)

uRange = obj.U;
x0Range = obj.R0;
u1=size(uRange,1); %shape of input ports
nInputs = size(allU, 2); % total number of inputs
if ~all(rad(x0Range) == 0) %find number of external inputs
    nExtInputs=nInputs-dim(x0Range);
else
    nExtInputs=nInputs;
end
assert(mod(nInputs, u1) == 0, 'total number of inputs is not divisible by number of input ports');
lowerBound=repelem(uRange.inf,round(nExtInputs/u1))';
upperBound = repelem(uRange.sup,round(nExtInputs/u1))';

if ~all(rad(x0Range) == 0) %append initial set bounds if not exact
    lowerBound=[lowerBound,x0Range.inf'];
    upperBound=[upperBound,x0Range.sup'];
end

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
MultiObj.nVar = nInputs;  % Set problem dimension
MultiObj.var_min = lowerBound;
MultiObj.var_max = upperBound;

%Fit Gaussian Process Meta Mod  el
GPmod = OK_Rmodel_kd_nugget(allU, allRob, 0, 2);

% optimize EI and CD with MOPSO
MultiObj.fun = @(x)[-EIcalc_kd(x,allU,GPmod,allRob), -CrowdingDist_kd(x,allU)];

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
newU = newSample(1:nExtInputs);
tak = (0:obj.ak.dt:obj.T)'; %define autokoopman time points
u= [tak,reshape(newU,[],u1)]; %append time points
x0 = newSample(nExtInputs+1:end);
end

