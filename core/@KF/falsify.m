function [solns,allDatas] = falsify(obj,varargin)

% falsify - Given a model and a set of specs (safe/unsafe set/stl),
%   perform the core falsification procedure to find a falsifying trajectory
%
% Syntax:
%    [solns, trainset] = falsify(obj)
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
%    solns - cell array of structs (equal in length to the number of runs)
%           containing the results of the falsification process, including
%           the falsifying trajectory, simulation time, and other information.
%
%    allData - Structure containing all data, including the trajectories and
%               inputs used during the falsification process as well as all
%               koopman models.
%
% See also:
%
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------
solns=cell(1,obj.runs);
allDatas=cell(1,obj.runs);
%initialize progress bar
if obj.verb==0
    msg = sprintf('KF runs completed: 0/%d',obj.runs);
    fprintf(msg);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
for run=1:obj.runs
    runtime=tic;
    %initialization and init assertions
    [obj,trainset,soln,specSolns,allData] = initialize(obj);
    trainIter = 0;
    falsified=false;
    perturb=0; %perturbation percentage for neighborhood reset
    tak = (0:obj.ak.dt:obj.T)'; %define autokoopman time points
    tsim = (0:obj.dt:obj.T)'; %time points for interpolating input
    %initialize simulation progress bar
    if obj.verb==1
        msg = sprintf('simulations exhausted: 0/%d',obj.maxSims);
        fprintf(msg);
        reverseSimStr = repmat(sprintf('\b'), 1, length(msg));
    end

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
            if (isnumeric(obj.nResets) && numel(trainset.X) >= obj.nResets) || (strcmp(obj.nResets,'auto') && trainIter>0 && getMinNormDistance(critX,critU,trainset,obj.R0,obj.U,obj.verb)<0.1)
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

            %if first iter, random or neighborhood training selected, or critical trajectory is
            %repeated, retrain with new xu else retrain with prev traj
            if trainIter==0 || obj.trainStrat>=1 || checkRepeatedTraj(critX,critU,trainset,obj.verb)
                if trainIter>0 && obj.solver.opts.usex0==1 && checkRepeatedTraj(critX,critU,trainset,obj.verb)
                    obj.solver.opts.usex0=0; %turn off warmstarting if repeated trajectory returned by solver
                    disp('Turned off warmstarting due to repeated solutions')
                end
                [t,x,u,simTime] = sampleSimulation(obj,allData,perturb);
                perturb=min(1,perturb+obj.sampPerturb); %increase perturbation
                soln.sims = soln.sims+1;
                soln.simTime = soln.simTime+simTime;
                %check if random input falsifies system, and break if it does
                [soln,falsified,robustness,~,newBest_]=checkFalsification(soln,x,u,t,obj.spec,obj.inputInterpolation,'reset simulation',obj.verb);
                allData.X{end+1}=x; allData.XU{end+1}=u; allData.t{end+1}=t; allData.Rob=[allData.Rob;robustness];
                if nargout>1;allData.koopModels{end+1}=[];end %store empty model as we are in reset
                if newBest_; perturb=obj.sampPerturb; end %reset pertrubation if new best soln found
                if falsified; break; end

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
        soln.koopTime = soln.koopTime+koopTime;
        % compute reachable set for Koopman linearized model (if reachability is used)
        if obj.reach.on
            reachTime=tic;
            try
                R = reachKoopman(obj,koopModel);
            catch
                vprintf(obj.verb,2,"error encountered whilst computing reachable set, resetting training data \n")
                trainIter=0;
                continue
            end
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
            vprintf(obj.verb,2,"error encountered whilst setup/solving, resetting training data \n")
            trainIter=0;
            continue;
        end

        %get critical spec with minimum robustness and corresponding soln struct
        [~,minIndex]=min(specSolns.values.rob);
        keys = specSolns.keys('cell');
        spec=keys{minIndex};
        curSoln=specSolns(spec);

        % this section check if critical trajectory is falsifying. If not, it also offsets if neccassary
        if curSoln.rob<inf %found a viable solution
            offsetIter = 0;
            while offsetIter <= max(obj.offsetStrat,0) %repeat this loop only if offset in same iteration is selected (offsetStrat=1)
                [critX0, critU] = falsifyingTrajectory(obj,curSoln);
                %interpolate input in accordance with interpolation strategy defined
                if ~isempty(critU)
                    assert(size(critU,1)>=2,'Input must have at least two sample points')
                    assert(size(critU,2)>=2,'Input must have at least two columns, where first column is time points')
                    usim = interp1(critU(:,1),critU(:,2:end),tsim,obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
                    usim =  max(obj.U.inf',min(obj.U.sup',usim)); %ensure that extrapolation is within input bounds
                    usim = [tsim,usim];
                else
                    usim=[]; %no input for the model
                end
                % run most critical inputs on the real system
                [t, critX, simTime] = simulate(obj, critX0, usim);
                soln.sims = soln.sims+1;
                soln.simTime = soln.simTime+simTime;

                %check if critical inputs falsify the system and store data
                [soln,falsified,robustness,Bdata,newBest_]=checkFalsification(soln,critX,critU,t,obj.spec,obj.inputInterpolation,'kf optimization',obj.verb);
                allData.X{end+1}=critX; allData.XU{end+1}=critU; allData.t{end+1}=t; allData.Rob=[allData.Rob;robustness];
                if nargout>1;allData.koopModels{end+1}=koopModel;end %store koop model if needed
                if newBest_; perturb=0; end %reset pertrubation if new best soln found
                if falsified; break; end

                if robustness~=inf %there exist a value for robustness for which we can neighborhood train or offset
                    if obj.trainStrat==1 && robustness >= soln.best.rob %neighborhood training mode
                        %remove last entry because it is not improving the obj
                        trainset.X(end) = []; trainset.XU(end)=[]; trainset.t(end) = [];
                    end
                    %if first offset iteration (re-solve with offset if offsetStrat==1 or save offset for next iter if offsetStrat==-1),
                    % robustness is greater than gap termination criteria for milp solver and an offset mode selected by user.
                    if offsetIter==0 && robustness > getMilpGap(obj.solver.opts) && abs(obj.offsetStrat)
                        assert(strcmp(spec.type,'logic'),'offset is currently only implemented for stl spec, please turn off offset by setting offsetStrat=0')
                        [critPreds]=bReachCulprit(Bdata,spec.set); %get predicates responsible for robustness value
                        if critPreds.Count > 0 %if there there exists predicates that are culprit for (+ve) robustness
                            obj.solver.opts.usex0=0; %avoid warmstarting if offsetting
                            Sys=specSolns(spec).koopMilp;
                            Sys.offsetMap = critPreds;
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
        %update simulations progress bar
        if obj.verb==1
            % Display the progress
            msg = sprintf('simulations exhausted: %d/%d',soln.sims,obj.maxSims); %Don't forget this semicolon
            fprintf([reverseSimStr, msg]);
            reverseSimStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    %close simulink obj
    close_system('ErrorIfShadowed',0);
    %assign solution result
    soln.falsified=falsified;
    soln.runtime=toc(runtime); %record runtime
    %remove simulations bar
    if obj.verb==1; fprintf(reverseSimStr); end

    LogicalStr = {'No', 'Yes'};
    vprintf(obj.verb,1,'<-------------------------------------------------------> \n')
    vprintf(obj.verb,1,"Run: %d/%d \n",run,obj.runs)
    vprintf(obj.verb,1,"Falsified: %s \n",LogicalStr{soln.falsified+1})
    vprintf(obj.verb,1,"number of simulations %d \n",soln.sims)
    vprintf(obj.verb,1,"Time taken %.2f seconds\n",soln.runtime)
    %store soln struct for this run
    solns{run}=soln;
    allDatas{run}=allData;
    %update progress bar
    if obj.verb==0
        % Display the progress
        msg = sprintf('KF runs completed: %d/%d',run,obj.runs); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
%remove progress bar
if obj.verb==0; fprintf(reverseStr); end
end

function repeatedTraj = checkRepeatedTraj(critX,critU,trainset,verb)
%check if critical initial set & input are the vsame as found before
repeatedTraj = false;
for r = 1:length(trainset.X)
    if isequal(critX(1,:)',trainset.X{r}(:,1)) && isequal(critU(:,2:end)',trainset.XU{r})
        repeatedTraj = true;
        vprintf(verb,2,2,"repeated critical trajectory, generating a new trajectory \n")
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
