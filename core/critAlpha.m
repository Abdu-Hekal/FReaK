function kfModel = critAlpha(R,kfModel)
% detemrine most critical reachable set and specification based on the
% robustnes

% loop over all specifications
spec=kfModel.spec;
rob = inf;

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
        %compute max time required to falsify stl and use it if less than
        %sim time (avoids unnecessary optim variables)
        maxStlSteps = min((maxStlTime(spec(i,1).set)/kfModel.ak.dt)+1,length(R.zono));
        %get prev solns
        prevSpecSol = kfModel.specSolns(kfModel.spec(i,1));
        %setup and run bluSTL
        tic
        %if no prev soln for this spec, setup alpha & stl milp vars/constrs
        try
            Sys=prevSpecSol.lti; %get previously setup milp problem with stl
            if ~kfModel.useOptimizer
                Sys=setupStl(Sys,true); %encode stl using milp
            end
        catch
            Sys=Koopman_lti(R.zono(1:maxStlSteps),kfModel.ak.dt,kfModel.solver.dt);
            if ~kfModel.pulseInput %if not pulse input, set cpBool
                Sys.cpBool=kfModel.cpBool;
            end
            Sys = setupAlpha(Sys);
            %convert disjunct stl from CORA format to blustl
            disjSet = disjunctiveNormalForm(spec(i,1).set); %is this necassary?
            bluStl = coraBlustlConvert(disjSet); %convert from cora syntax to blustl
            Sys.stlList = {bluStl};
            if kfModel.useOptimizer
                Sys=setupStl(Sys,false); %encode stl using milp
            else
                Sys=setupStl(Sys,true);
            end
        end

        Sys.reachZonos=R.zono(1:maxStlSteps); %update reach zonos with new
        Sys = setupReach(Sys);
        kfModel.soln.milpSetupTime = kfModel.soln.milpSetupTime+toc;

        tic
        if kfModel.useOptimizer
            Sys = setupOptimizer(Sys,kfModel.solver.opts);
        end
        specSoln.lti=Sys; %store lti object with setcup optimizer
        Sys=optimize(Sys,kfModel.solver.opts);
        kfModel.soln.milpSolvTime =kfModel.soln.milpSolvTime+toc;

        %get results
        rob_ = value(Sys.Pstl);
        alpha = value(Sys.Alpha);

        %clear solution from yalmip and assign only alpha for warmstarting
        if kfModel.solver.opts.usex0
            yalmip('clearsolution')
        end

        %TODO: how can we compare stl robustness and reachset robustness.
        if rob_ < rob
            alphaCrit = alpha';
            setCrit = R.set{end};
            rob = rob_;
            specCrit=spec(i,1);
        end
    else
        error('This type of specification is not supported!');
    end

    %store solution for this iteration for each spec.
    specSoln.rob=rob_; specSoln.alpha=alpha;
    kfModel.specSolns(kfModel.spec(i,1)) = specSoln;
end
%store solution for this iteration
kfModel.soln.rob=rob; kfModel.soln.alpha=alphaCrit;
kfModel.soln.set=setCrit; kfModel.soln.spec=specCrit;
end

function [r,alpha] = robustness(P,Z)
% compute robustness of the zonotope Z with respect to an unsafe polytope P

% catch special case of a halfspace to accelearte computation
%     if size(P.P.A,1) == 1
%         disp("halfspace")
%     end
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
    function maxTime = maxStlTime(stl)
        maxTime=0;
        timedOp = {'finally', 'globally', 'release', 'until'};
        isMember = ismember(stl.type, timedOp);

        if isMember
            maxTime=maxTime+obj.to;
        end
        if ~isempty(obj.lhs)
            maxTime=maxTime(obj.lhs);
        end
        if ~isempty(obj.rhs)
            maxTime=maxTime(obj.rhs);
        end
    end

end
