function [soln,trainset] = falsify(obj,varargin)

% falsify - Given a model and a set of specs (safe/unsafe set/stl),
%   perform the core falsification procedure to find a falsifying trajectory
%
% Syntax:
%    [obj, trainset] = falsify(obj)
%
% Description:
%    This function is the core of the falsification procedure, which iteratively
%    searches for a falsifying trajectory by simulating the system, training
%    the Koopman model, and optimizing for critical trajectories. The process
%    continues until a falsifying trajectory is found or a specified timeout is reached.
%
% Inputs:
%    obj - KF object containing the Koopman model and various parameters
%              needed for the falsification process.
%    trainset (optional) - initial trainset to warmstart training
%
% Outputs:
%    obj - KF object containing the results of the falsification process,
%              including the falsifying trajectory, simulation time, and other information.
%
%    trainset - Structure containing the training data, including the trajectories and
%               inputs used during the falsification process.
%
% See also:
%
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


runtime=tic;
%initialization and init assertions
[obj,trainset,soln,specSolns] = initialize(obj);
trainIter = 0;
falsified=false;

if nargin > 1
    trainset = varargin{1};
    validateTrainset(trainset,obj.U)
end
while soln.sims <= obj.maxSims && ~falsified
    if ~(soln.sims==0 && ~isempty(trainset.t)) %jump straight to learning if an initial trainset is provided
        %timeout
        if toc(runtime) > obj.timeout
            break
        end
        %reset after size of trainset==nResets;
        if (isnumeric(obj.nResets) && numel(trainset.X) >= obj.nResets) || (strcmp(obj.nResets,'auto') && trainIter>0 && getMinNormDistance(critX,critU,trainset,obj.R0,obj.U,obj.verb)<0.01)
            trainIter = 0;
            %reset offsets
            for ii=1:numel(obj.spec)
                spec=obj.spec(ii);
                specSolns(spec).koopMilp.offsetMap=containers.Map('KeyType', 'double', 'ValueType', 'double');
            end
            vprintf(obj.verb,2,2,'Reset applied to training set of size %d\n',size(trainset.t,2))
        end
        % empty trainset at reset, or empty trainset after first iter because it is random trajectory, if rmRand is selected by user.
        if  trainIter ==0 || (trainIter == 1 && obj.rmRand)
            trainset.X = {}; trainset.XU={}; trainset.t = {};
        end

        %if first iter, random trajectory setting selected, or critical trajectory is
        %repeated, retrain with random xu else retrain with prev traj
        if trainIter==0 || obj.trainRand>=2 || ( obj.trainRand==1 && rem(trainIter, 2) == 0) || checkRepeatedTraj(critX,critU,trainset,obj.verb)
            if trainIter>0 && obj.solver.opts.usex0==1 && checkRepeatedTraj(critX,critU,trainset,obj.verb)
                obj.solver.opts.usex0=0; %turn off warmstarting if repeated trajectory returned by solver
                disp('Turned off warmstarting due to repeated solutions')
            end
            [t, x, u, simTime] = randSimulation(obj);
            soln.sims = soln.sims+1;
            soln.simTime = soln.simTime+simTime;
            %check if random input falsifies system, and break if it does
            [soln,falsified]=checkFalsification(soln,x,u,t,obj.spec,obj.inputInterpolation,'reset simulation',obj.verb);
            if falsified; break; end

            tak = (0:obj.ak.dt:obj.T)'; %define autokoopman time points
            xak = interp1(t,x,tak,obj.trajInterpolation); %define autokoopman trajectory points
        else
            xak=interp1(t,critX,tak,obj.trajInterpolation); %pass x0 as full x to avoid simulation again
            u=critU;
        end

        %add trajectory to koopman trainset
        trainset.t{end+1} = tak;
        trainset.X{end+1} = xak';
        trainset.XU{end+1} = u(:,2:end)';

    end

    %run autokoopman and learn linearized model
    [koopModel,koopTime] = learnKoopModel(obj, trainset);
    soln.koopModel=koopModel; %store koopman model
    soln.koopTime = soln.koopTime+koopTime;
    % compute reachable set for Koopman linearized model (if reachability is used)
    if obj.reach.on
        reachTime=tic;
        R = reachKoopman(obj,koopModel);
        soln.reachTime=soln.reachTime+toc(reachTime); %time for reachability computation
    else
        R=[];
    end
    % determine most critical reachable set and specification
    try
        optimTime=tic;
        specSolns = critAlpha(obj,R,koopModel,specSolns);
        soln.optimTime=soln.optimTime+toc(optimTime);
    catch
        disp("error encountered whilst setup/solving, resetting training data")
        trainIter=0;
        continue;
    end

    %get critical spec with minimum robustness and corresponding soln struct
    [~,minIndex]=min(specSolns.values.rob);
    keys = specSolns.keys('cell');
    spec=keys{minIndex};
    curSoln=specSolns(spec);

    if curSoln.rob<inf %found a viable solution
        offsetIter = 0;
        while offsetIter <= max(obj.offsetStrat,0) %offset once if offset in same iteration is selected (offsetStrat=1)
            [critX0, critU] = falsifyingTrajectory(obj,curSoln);

            % run most critical inputs on the real system
            [t, critX, simTime] = simulate(obj, critX0, critU);
            soln.sims = soln.sims+1;
            soln.simTime = soln.simTime+simTime;

            [soln,falsified,robustness,Bdata]=checkFalsification(soln,critX,critU,t,obj.spec,obj.inputInterpolation,'kf optimization',obj.verb);
            if falsified; break; end

            if robustness~=inf %there exist a value for robustness for which we can neighborhood train or offset
                if obj.trainRand==2 && robustness >= soln.best.rob %neighborhood training mode
                    %remove last entry because it is not improving the obj
                    trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
                end
                %if first offset iteration (re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1),
                % robustness is greater than gap termination criteria for milp solver and an offset mode selected by user.
                if offsetIter==0 && robustness > getMilpGap(obj.solver.opts) && abs(obj.offsetStrat)
                    assert(strcmp(spec.type,'logic'),'offset is currently only implemented for stl spec, please turn off offset by setting offsetStrat=0')
                    offsetMap=bReachCulprit(Bdata,spec.set); %get predicates responsible for robustness value
                    if offsetMap.Count > 0 %if there there exists predicates that are culprit for (+ve) robustness
                        obj.solver.opts.usex0=0; %avoid warmstarting if offsetting
                        Sys=specSolns(spec).koopMilp;
                        Sys.offsetMap = offsetMap;
                        if obj.offsetStrat == 1 %if offset strategy in this iteration selected
                            set = spec.set;
                            if ~isequal(set,Sys.stl) || ~obj.solver.useOptimizer
                                Sys.stl = set;
                                Sys=setupStl(Sys,~obj.solver.useOptimizer); %encode stl using milp
                            end
                            Sys=optimize(Sys,obj.solver.opts);
                            curSoln.alpha = value(Sys.alpha); %new alpha value after offset
                            curSoln.u = value(Sys.u);
                            % TODO: if offset gives better val of robustness, should we pass
                            %                     % it as training data instead? should we pass both?
                        else %obj.offsetStrat == -1: offset next iteration
                            specSolns(spec).koopMilp=Sys;
                        end
                    else
                        break
                    end
                else
                    break
                end
            end
            offsetIter = offsetIter+1;
        end
    end
    %clear previous solution from yalmip (makes warmstart feasible)
    if obj.solver.opts.usex0
        yalmip('clearsolution')
    end
    trainIter=trainIter+1;
end
%close simulink obj
close_system('ErrorIfShadowed',0);
%assign solution result
soln.falsified=falsified;
soln.runtime=toc(runtime); %record runtime

LogicalStr = {'No', 'Yes'};
vprintf(obj.verb,1,'<-------------------------------------------------------> \n')
vprintf(obj.verb,1,"Falsified: %s \n",LogicalStr{soln.falsified+1})
vprintf(obj.verb,1,"number of simulations %d \n",soln.sims)
vprintf(obj.verb,1,"Time taken %.2f seconds\n",soln.runtime)
end

function repeatedTraj = checkRepeatedTraj(critX,critU,trainset,verb)
%check if critical initial set & input are the same as found before
repeatedTraj = false;
for r = 1:length(trainset.X)

    %     A=trainset.XU{r};
    %     B=critU(:,2:end)';
    %     normalized_distance = norm(A - B) / (sqrt(2) * max(norm(A), norm(B)));
    %     vprintf(verb,2,"normalize distance with trainset %d is %f \n",r,normalized_distance)

    if isequal(critX(1,:)',trainset.X{r}(:,1)) && isequal(critU(:,2:end)',trainset.XU{r})
        repeatedTraj = true;
        vprintf(verb,3,2,"repeated critical trajectory, generating a new random trajectory \n")
        break
    end
end
end

function minND=getMinNormDistance(critX,critU,trainset,R0,U,verb)
minND=inf;
for r = 1:length(trainset.X)

    A=[];
    B=[];
    if ~isempty(R0) && any(rad(R0) ~= 0) %non-exact initial set
        range=R0.sup-R0.inf;
        nonzero=find(range);
        A = (trainset.X{r}(nonzero,1)-R0.inf(nonzero))./(range(nonzero));
        B = (critX(1,nonzero)'-R0.inf(nonzero))./(range(nonzero));
    end
    if ~isempty(U) && any(rad(U) ~= 0) %non-exact inputs
        range=U.sup-U.inf;
        nonzero=find(range);
        a=reshape(trainset.XU{r}(nonzero,:)./(range(nonzero)),[],1);
        A = [A;a];
        b=critU(:,2:end)';
        b=b(nonzero,:);
        b=reshape(b./(range(nonzero)),[],1);
        B=[B;b];
    end
    ND = norm(A - B)/numel(A);
    if ND < minND
        minND=ND;
    end
end
vprintf(verb,3,"min normalize distance with trainset is %f \n",minND)
end

function gap=getMilpGap(opts)
% Note that if robustness is less than gap, offset most likely is not benefecial.
if opts.solver == "gurobi"
    gap=opts.gurobi.MIPGapAbs;
elseif opts.solver == "cplex"
    gap=opts.cplex.epagap;
else
    gap=0.1;
end
end

function validateTrainset(trainset,U)
assert(isstruct(trainset), 'trainset must be a struct');
assert(isfield(trainset, 'X'), 'trainset must have a field named "X" which is a cell of training trajectories');
assert(isfield(trainset, 't'), 'trainset must have a field named "t" which is a cell of times corresponding to training trajectories');
assert(numel(trainset.X)==numel(trainset.t),'number of training trajectories must be the same as number of time arrays')
for ii=1:numel(trainset.t)
    assert(size(trainset.X{ii},2)==size(trainset.t{ii},1),'length of training trajectory must be same as time points')
end
if ~isempty(U)
    assert(isfield(trainset, 'XU'), 'model has inputs, trainset must have a field named "XU" which is a cell of training inputs');
    assert(numel(trainset.XU)==numel(trainset.t),'number of input training data must be the same as number of time arrays and training trajectories')
    for ii=1:numel(trainset.t)
        assert(size(trainset.XU{ii},2)==size(trainset.t{ii},1),'length of training inputs must be same as time points')
    end
end
end
