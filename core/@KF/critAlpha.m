function specSolns = critAlpha(obj,R,koopModel,specSolns)
% critAlpha - Determine the most critical alpha values (or inputs u) based
% on robustness. This is usually achieved by solving an optimization unless
% the unsafe/safe set is a halfspace.
%
% Syntax:
%    specSolns = critAlpha(obj,R, koopModel)
%
% Description:
%    This function iterates over all specifications and computes the most
%    critical alpha values (or inputs) based on robustness measures.
%    The result is stored in the obj object, updating the critical
%    alpha values, control input, reachable set, and specification.
%
% Inputs:
%    obj - KF object containing the Koopman model, specifications, and
%              various parameters needed for the falsification process.
%    R - Reachable set information, including sets, time intervals, and
%        zonotopes.
%    koopModel - struct representing koopman model with following:
%       A - State transition matrix of the Koopman linearized model.
%       B - Input matrix of the Koopman linearized model.
%       g - observables function of the Koopman linearized model.
%   specSolns - soln dictionary for each spec from previous iteration
%
% Outputs:
%    specSolns - soln dictionary for each spec with critical information, including
%              alpha values, control input, reachable set.
%
%
% See also: falsify
%
% Author:      Niklas Kochdumper, Abdelrahman Hekal
% Written:     28-February-2023
% Last update: 4-December-2023
% Last revision: [Date]
%
% -------------------------- Auxiliary Functions --------------------------
%
% robustness - Compute robustness of a zonotope with respect to an unsafe
%              polytope.
%
% Syntax:
%    [r, alpha] = robustness(P, Z)
%
% Description:
%    This function computes the robustness of a zonotope Z with respect to
%    an unsafe polytope P. It is used as an auxiliary function in critAlpha.
%
% Inputs:
%    P - Unsafe polytope information.
%    Z - Zonotope information.
%
% Outputs:
%    r - Robustness measure.
%    alpha - Critical alpha values.
%
% See also: critAlpha
%
% Author:      Niklas Kochdumper, Abdelrahman Hekal
% Written:     28-February-2023
% Last update: [Date]
% Last revision: [Date]
%
% -----------------------------------------------------------------------------
% ------------- BEGIN CODE --------------

% loop over all specifications
spec=obj.spec;
for i = 1:size(spec,1)
    rob = inf; %initial robustness value
    setCrit=[];
    u=[]; %initialize empty list to store critical inputs
    x=[]; %initialize empty array to store crit traj
    specSoln = struct; %struct for stored solution

    % different types of specifications
    if strcmp(spec(i,1).type,'unsafeSet')

        set = mptPolytope(spec(i,1).set);

        % loop over all reachable sets
        for j = 1:length(R.set)
            if isempty(spec(i,1).time) || ...
                    isIntersecting(R.time{i},spec(i,1).time)

                % compute robustness
                [rob_,alpha_] = robustness(set,R.zono{j});

                % check if smaller than current robustness
                if rob_ < rob
                    alpha=alpha_;
                    setCrit = R.set{j};
                    rob = rob_;
                end
            end
        end

    elseif strcmp(spec(i,1).type,'safeSet')

        set = mptPolytope(spec(i,1).set);

        % loop over all reachable sets
        for j = 1:length(R.set)
            if isempty(spec(i,1).time) || ...
                    isIntersecting(R.time{i},spec(i,1).time)

                % compute robustness
                dist = infimum(interval(-set.P.A*R.zono{j})) + set.P.b;
                [rob_,ind] = max(dist);
                alpha_ = -sign(-set.P.A(ind,:)*generators(R.zono{j}));

                % check if smaller than current robustness
                if rob_ < rob
                    alpha=alpha_;
                    setCrit = R.set{j};
                    rob = rob_;
                end
            end
        end

    elseif strcmp(spec(i,1).type,'logic')
        %get prev solns
        prevSpecSol = specSolns(obj.spec(i,1));
        %setup and run solver
        tic
        %if no prev soln for this spec, setup alpha & stl vars/constrs
        if isfield(prevSpecSol, 'KoopSolver')
            Sys=prevSpecSol.KoopSolver; %get previously setup solver with stl
        else
            Sys=KoopSolver(obj.T,obj.ak.dt,obj.solver.timePoints,obj.R0,obj.U);
            Sys.normalize = obj.solver.normalize; %set normalization setting
            if ~obj.pulseInput %if not pulse input, set cpBool
                Sys.cpBool=obj.cpBool;
            end
            if obj.reach.on
                Sys.reachZonos=R.zono; %update reach zonos with new
                Sys = setupAlpha(Sys);
            else
                Sys.nObs = obj.ak.nObs;
                Sys = setupInit(Sys);
            end
            Sys = setupCP(Sys); %setup control points constraints
        end

        %setup problem from scratch if number of generators is no longer
        %the same.
        if obj.reach.on
            if size(Sys.alpha,2) ~= size(generators(R.zono{end}),2)
                Sys=KoopSolver(obj.T,obj.ak.dt,obj.solver.timePoints,obj.R0,obj.U);
                Sys.normalize = obj.solver.normalize; %set normalization setting
                if ~obj.pulseInput %if not pulse input, set cpBool
                    Sys.cpBool=obj.cpBool;
                end
                Sys.reachZonos=R.zono; %update reach zonos with new
                Sys = setupAlpha(Sys);
            end
        end

        %if auto contraints, add constraints for critical times and
        %corresponding predicates, also check if critTimes have been
        %defined
        set=spec(i,1).set;
        if obj.solver.autoAddConstraints==2
            %initialize weights if they have not yet been set
            if isempty(Sys.weights)
                %first get all predicstes
                breachPhi = STL_Formula('phi',coraBreachConvert(set));
                mus = STL_ExtractPredicates(breachPhi);
                %create a weight for every step in ak
                Sys.weights=ones(numel(mus),(obj.T/obj.ak.dt)+1);
            end
            if ~isempty(prevSpecSol.critTimes)
                %add weights for new constrs
                for c=1:numel(prevSpecSol.critTimes)
                    constr=prevSpecSol.critTimes{c};
                    p = find(strcmp(prevSpecSol.preds, constr.pred), 1);
                    %we increase the weight of each critical time by how
                    %much it violates the stl spec
                    Sys.weights(p,constr.time+1)=Sys.weights(p,constr.time+1)+constr.weight; %+constr.weight
                end
            end
            %setup weighted stl from scratch: if (1) first time encoding stl OR
            % (2) weight is hardcoded every time OR
            % (3) time points for evaluating stl changes, typically for 'auto' solver step
            if ~isequal(set,Sys.stl) || ~obj.solver.useOptimizer || ~isequal(obj.solver.timePoints,Sys.solverTimePoints)
                Sys.stl = set;
                Sys.solverTimePoints=obj.solver.timePoints;
                Sys=setupStl(Sys,~obj.solver.useOptimizer,true); %encode stl using weighted pred constraints
            end
            specSoln.critTimes={}; %empty list for crit times to only store new crit Times
            specSoln.preds=prevSpecSol.preds;
        elseif obj.solver.autoAddConstraints==1
            Sys.solverTimePoints=obj.solver.timePoints;
            if ~isempty(prevSpecSol.critTimes)
                Sys=addPredConstr(Sys,prevSpecSol.critTimes,prevSpecSol.preds,~obj.solver.useOptimizer,obj.offsetStrat);
            end
            specSoln.critTimes=prevSpecSol.critTimes; %store all crit Times
            specSoln.preds=prevSpecSol.preds;
        elseif obj.solver.autoAddConstraints==0
            %setup stl from scratch: if (1) first time encoding stl OR
            % (2) offset strategy is used and optimizer object is not used, i.e. offset is hardcoded
            % every time OR (3) time points for evaluating stl changes, typically for 'auto' solver step
            if ~isequal(set,Sys.stl) || (~obj.solver.useOptimizer && numEntries(Sys.offsetMap)>0) || ~isequal(obj.solver.timePoints,Sys.solverTimePoints)
                Sys.stl = set;
                Sys.solverTimePoints=obj.solver.timePoints;
                Sys=setupStl(Sys,~obj.solver.useOptimizer,false); %encode stl using milp
            end
        end

        %setup evolution of system using reachable sets or direct encoding
        if obj.reach.on
            Sys.reachZonos=R.zono; %update reach z onos with new
            Sys = setupReach(Sys);
        else
            Sys.A=koopModel.A;
            Sys.B=koopModel.B;
            Sys.g=koopModel.g;
            Sys=setupDynamics(Sys);
        end

        %setup optimizer object
        if obj.solver.useOptimizer
            Sys = setupOptimizer(Sys,obj.solver.opts);
        end
        specSoln.KoopSolver=Sys; %store KoopSolver object with setup optimizer
        Sys=optimize(Sys,obj.solver.opts);

        %get results
        rob = value(Sys.Pstl);
        alpha = value(Sys.alpha);
        u = value(Sys.u);
        x=value(Sys.x);

        if ~isempty(R)
            setCrit = R.set{end};
        end
    else
        error('This type of specification is not supported!');
    end

    %TODO: how can we compare stl robustness and reachset robustness.
    %store solution for this iteration for each spec.
    specSoln.rob=rob; specSoln.alpha=alpha; specSoln.u=u; specSoln.x=x;
    specSoln.set=setCrit;
    specSolns(spec(i,1)) = specSoln;
end
end

% -------------------------- Auxiliary Functions --------------------------

function [r,alpha] = robustness(P,Z)
% compute robustness of the zonotope Z with respect to an unsafe polytope P

% catch special case of a halfspace to accelerate computation
if size(P.P.A,1) == 1
    r = infimum(interval(P.P.A*Z)) - P.P.b;
    alpha = -sign(P.P.A*generators(Z))';
else
    if isIntersecting(P,Z)

        % solve linear program with variables [x;r;\alpha]
        %
        %   max r s.t. A(i,:)*x + ||A(i,:)||*r <= b,
        %              x = c + G * \alpha,
        %              r >= 0,
        %              -1 <= \alpha <= 1

        A = P.P.A; b = P.P.b; n = size(A,2);
        c = center(Z); G = generators(Z); m = size(G,2);

        % constraint A(i,:)*x+||A(i,:)||*r <= b
        C1 = [A,sum(A.^2,2)]; d1 = b;

        % constraint r >= 0
        C2 = [zeros(1,size(A,2)),-1]; d2 = 0;

        % constraint -1 <= \alpha <= 1
        C3 = [eye(m);-eye(m)]; d3 = ones(2*m,1);

        % constraint x = c + G * \alpha,
        Ceq = [eye(n),zeros(n,1),-G]; deq = c;

        % combined inequality constraints
        C = blkdiag([C1;C2],C3); d = [d1;d2;d3];

        % objective function
        f = zeros(size(Ceq,2),1); f(n+1) = -1;

        % solve linear program
        options = optimoptions('linprog','display','off');

        [x,r] = linprog(f,C,d,Ceq,deq,[],[],options);

        alpha = x(n+2:end);

    else

        % solve quadratic program with variables [d;x;\alpha]
        %
        %   min ||d||^2 s.t. d = c + G*\alpha - x,
        %                    A*x <= b,
        %                    -1 <= \alpha <= 1

        A = P.P.A; b = P.P.b; n = size(A,2);
        c = center(Z); G = generators(Z); m = size(G,2);

        % constraint d = c + G*\alpha - x
        Ceq = [eye(n),eye(n),-G]; deq = c;

        % constraint A*x <= b
        C1 = [zeros(size(A,1),n),A]; d1 = b;

        % constraint -1 <= \alpha <= 1
        C2 = [eye(m);-eye(m)]; d2 = ones(2*m,1);

        % combined inequality constraints
        C = blkdiag(C1,C2); d = [d1;d2];

        % objective function
        H = blkdiag(2*eye(n),zeros(n+m));
        f = zeros(2*n+m,1);

        % solve quadratic program
        options = optimoptions('quadprog','display','off');

        [x,r] = quadprog(H,f,C,d,Ceq,deq,[],[],[],options);

        r = sqrt(r);
        alpha = x(2*n+1:end);
    end
end
end
