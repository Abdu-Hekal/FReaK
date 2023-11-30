function obj = critAlpha(obj,R,koopModel)
% critAlpha - Determine the most critical alpha values (or inputs u) based 
% on robustness. This is usually achieved by solving an optimization unless
% the unsafe/safe set is a halfspace. 
%
% Syntax:
%    obj = critAlpha(obj,R, koopModel)
%
% Description:
%    This function iterates over all specifications and computes the most
%    critical alpha values (or inputs) based on robustness measures.
%    The result is stored in the obj object, updating the critical
%    alpha values, control input, reachable set, and specification.
%
% Inputs:
%    R - Reachable set information, including sets, time intervals, and
%        zonotopes.
%    koopModel - struct representing koopman model with following:
%       A - State transition matrix of the Koopman linearized model.
%       B - Input matrix of the Koopman linearized model.
%       g - observables function of the Koopman linearized model.
%    obj - KF object containing the Koopman model, specifications, and
%              various parameters needed for the falsification process.
%
% Outputs:
%    obj - Updated KF object with critical information, including
%              alpha values, control input, reachable set, and specification.
%
%
% See also: falsify
%
% Author:      Niklas Kochdumper, Abdelrahman Hekal
% Written:     28-February-2023
% Last update: [Date]
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
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: [Date]
% Last revision: [Date]
%
% -----------------------------------------------------------------------------
% ------------- BEGIN CODE --------------

% loop over all specifications
spec=obj.spec;
rob = inf;
uCrit = []; %initialise critical input.

for i = 1:size(spec,1)

    %struct for stored solution
    specSoln = struct;

    % different types of specifications
    if strcmp(spec(i,1).type,'unsafeSet')

        set = mptPolytope(spec(i,1).set);

        % loop over all reachable sets
        for j = 1:length(R.set)
            if isempty(spec(i,1).time) || ...
                    isIntersecting(R.time{i},spec(i,1).time)

                % compute robustness
                [rob_,alpha] = robustness(set,R.zono{j});

                % check if smaller than current robustness
                if rob_ < rob
                    alphaCrit = alpha;
                    setCrit = R.set{j};
                    rob = rob_;
                    specCrit=spec(i,1);
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
                alpha = -sign(-set.P.A(ind,:)*generators(R.zono{j}));

                % check if smaller than current robustness
                if rob_ < rob
                    alphaCrit = alpha;
                    setCrit = R.set{j};
                    rob = rob_;
                    specCrit=spec(i,1);
                end
            end
        end

    elseif strcmp(spec(i,1).type,'logic')
        %get prev solns
        prevSpecSol = obj.specSolns(obj.spec(i,1));
        %setup and run milp
        tic
        %if no prev soln for this spec, setup alpha & stl milp vars/constrs
        try
            Sys=prevSpecSol.lti; %get previously setup milp problem with stl
        catch
            Sys=KoopMILP(obj.T,obj.ak.dt,obj.solver.dt,obj.R0,obj.U);
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
        %the same. TODO: check and clean this section
        if obj.reach.on
            if size(Sys.alpha,2) ~= size(generators(R.zono{end}),2)
                Sys=KoopMILP(obj.T,obj.ak.dt,obj.solver.dt,obj.R0,obj.U);
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
            end
        end

        set=spec(i,1).set;
        if ~isequal(set,Sys.stl) || (~obj.solver.useOptimizer && Sys.offsetMap.Count>0)
            Sys.stl = set;
            Sys=setupStl(Sys,~obj.solver.useOptimizer); %encode stl using milp
        end

        if obj.reach.on
            Sys.reachZonos=R.zono; %update reach zonos with new
            Sys = setupReach(Sys);
        else
            Sys.A=koopModel.A;
            Sys.B=koopModel.B;
            Sys.g=koopModel.g;
            Sys=setupDynamics(Sys);
        end

        obj.soln.milpSetupTime = obj.soln.milpSetupTime+toc;

        tic
        if obj.solver.useOptimizer
            Sys = setupOptimizer(Sys,obj.solver.opts);
        end
        specSoln.lti=Sys; %store lti object with setcup optimizer
        Sys=optimize(Sys,obj.solver.opts);
        obj.soln.milpSolvTime =obj.soln.milpSolvTime+toc;

        %get results
        rob_ = value(Sys.Pstl);
        alpha = value(Sys.alpha);
        u = value(Sys.u);


        %TODO: how can we compare stl robustness and reachset robustness.
        if rob_ < rob
            alphaCrit = alpha';
            uCrit = u;
            if ~isempty(R)
                setCrit = R.set{end};
            else
                setCrit=[];
            end
            rob = rob_;
            specCrit=spec(i,1);
        end
    else
        error('This type of specification is not supported!');
    end

    %store solution for this iteration for each spec.
    specSoln.rob=rob_; specSoln.alpha=alpha;
    obj.specSolns(obj.spec(i,1)) = specSoln;
end
%store solution for this iteration
obj.soln.rob=rob; obj.soln.alpha=alphaCrit; obj.soln.u=uCrit;
obj.soln.set=setCrit; obj.soln.spec=specCrit;
end

% -------------------------- Auxiliary Functions --------------------------

function [r,alpha] = robustness(P,Z)
% compute robustness of the zonotope Z with respect to an unsafe polytope P

% catch special case of a halfspace to accelearte computation
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
