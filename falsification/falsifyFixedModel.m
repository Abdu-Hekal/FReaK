function [x0,u] = falsifyFixedModel(A,B,g,dt,spec,R0,U,tFinal)
% determine the most critical initial point and input trajectory of a
% Koopman linearized model for the given trajectory
%
% Input arguments:
%
%   -A:        system matrix for the Koopman linearized model
%   -B:        input matrix for the Koopman linearized model
%   -g:        function handle to the observable function for the Koopman
%              linearized model
%   -dt:       time step size for the Koopman lineaized model
%   -spec:     specification defined as an object of the CORA specification
%              class (includes safe sets, unsafe sets, and temporal logic)
%   -R0:       initial set (CORA class interval)
%   -U:        input set (CORA class interval or zonotope)
%   -tFinal:   final time
%
% Output arguments:
%   
%   -x0:       initial state for the most critical trajectory
%   -u:        piecewise constant inputs for the most critical trajectory

    % compute reachable set for Koopman linearized model
    R = reachKoopman(A,B,g,R0,U,tFinal,dt);

    % determine most critical reachable set and specification
    [set,alpha] = mostCriticalReachSet(R,spec);

    % extract most critical initial state and input signal
    [x0,u] = falsifyingTrajectory(R0,U,set,alpha);

end


% Auxiliary Functions -----------------------------------------------------

function R = reachKoopman(A,B,g,R0,U,tFinal,dt)
% compute reachable set for Koopman linearized model

    % compute initial set using Taylor model arithmetic
    n = dim(R0);
    tay = taylm(R0);
    tay = g(tay);
    R0 = polyZonotope(tay);
    R0 = polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:));

    % compute reachable set
    t = 0:dt:ceil(tFinal/dt)*dt; U = zonotope(U);

    set = cell(length(t),1); time = set; zono = set;

    set{1} = R0; time{1} = interval(-dt/2,dt/2);

    for i = 1:length(t)-1
        set{i+1} = A*set{i} + B*zonotope(U);
        time{i+1} = time{i} + dt;
    end

    % compute ouput set
    for i = 1:length(set)
        set{i} = project(set{i},1:n);
        zono{i} = zonotope(set{i});
    end

    R.set = set; R.time = time; R.zono = zono;
end

function [setCrit,alphaCrit] = mostCriticalReachSet(R,spec)
% detemrine most critical reachable set and specification based on the
% robustnes

    % loop over all specifications
    rob = inf; 

    for i = 1:size(spec,1)
        
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
                    end
                end
            end

        elseif strcmp(spec(i,1).type,'logic')
    
            % convert to negation normal form
            phi = negationNormalForm(spec(i,1).set);

            % constrct time vector
            time = zeros(size(R.time));

            for j = 1:length(R.time)
                time(j) = center(R.time{j});
            end

            % compute robustness with recursive function
            [rob_,alpha,ind] = robustnessTemporalLogic(phi,R.zono,time);

            % check if smaller than current robustness
            if rob_(1) < rob
                alphaCrit = alpha{1};
                setCrit = R.set{ind(1)};
                rob = rob_(1);
            end

        else
            error('This type of specification is not supported!');
        end
    end
end

function [x0,u] = falsifyingTrajectory(R0,U,set,alpha)
% extract the most critical initial state and input signal from the most
% critical reachable set and specification

    % determine most critical initial state
    alphaInit = zeros(size(set.expMat,1),1);

    temp = prod(ones(size(set.expMat))-mod(set.expMat,2),1);
    expMat = set.expMat(:,temp == 0);

    ind = find(sum(expMat,1) == 1);

    for i = 1:size(set.expMat,1)
        for j = ind
            if expMat(i,j) == 1
                alphaInit(i) = alpha(j);
            end
        end
    end

    R0 = zonotope(R0);
    x0 = center(R0) + generators(R0)*alphaInit;

    % determine most ctritical control input
    if ~isempty(set.Grest)
    
        alphaU = alpha(size(set.G,2)+1:end);

        U = zonotope(U); c_u = center(U); G_u = generators(U);
        alphaU = reshape(alphaU,[size(G_u,2),length(alphaU)/size(G_u,2)]);

        u = c_u + G_u*alphaU;
    else
        u = center(U);
    end
end

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

function [r,alpha,ind] = robustnessTemporalLogic(phi,R,time)
% recursive function to compute the robustness of a temporal logic formula

    if ~phi.temporal

        % convert logic equation to union of safe sets
        eq = disjunctiveNormalForm(phi);
        safeSet = getClauses(eq,'dnf');

        for k = 1:length(safeSet)
            safeSet{k} = convert2set(safeSet{k});
        end

        % convert to a union of unsafe sets
        unsafeSet = safe2unsafe(safeSet);

        % loop over all reachable sets
        r = inf*ones(1,length(R)); alpha = cell(1,length(R));
        ind = 1:length(R);

        for i = 1:length(R)
            
            % loop over all unsafe sets
            for j = 1:length(unsafeSet)
                
                % compute robustness
                [r_,alpha_] = robustness(unsafeSet{j},R{i});

                % check if robustness is smaller than for other sets
                if r_ < r(i)
                    r(i) = r_; alpha{i} = alpha_;
                end
            end
        end

    elseif strcmp(phi.type,'&')

        [r1,alpha1,ind1] = robustnessTemporalLogic(phi.lhs,R,time);
        [r2,alpha2,ind2] = robustnessTemporalLogic(phi.rhs,R,time);

        r = zeros(size(r1)); ind = zeros(size(ind1)); 
        alpha = cell(size(alpha1)); 

        for i = 1:length(r)
            if r1(i) < r2(i)
                r(i) = r1(i); ind = ind1(i); alpha{i} = alpha1{i};
            else
                r(i) = r2(i); ind = ind2(i); alpha{i} = alpha2{i};
            end
        end

    elseif strcmp(phi.type,'|')

        [r1,alpha1,ind1] = robustnessTemporalLogic(phi.lhs,R,time);
        [r2,alpha2,ind2] = robustnessTemporalLogic(phi.rhs,R,time);

        r = zeros(size(r1)); ind = zeros(size(ind1)); 
        alpha = cell(size(alpha1)); 

        for i = 1:length(r)
            if r1(i) > r2(i)
                r(i) = r1(i); ind = ind1(i); alpha{i} = alpha1{i};
            else
                r(i) = r2(i); ind = ind2(i); alpha{i} = alpha2{i};
            end
        end

    elseif strcmp(phi.type,'next')

        [r,alpha,ind] = robustnessTemporalLogic(phi.lhs,R,time);
        index = find(time >= phi.from);
        r = r(index); alpha = alpha(index); ind = ind(index);

    elseif strcmp(phi.type,'finally')

        [r_,alpha_,ind_] = robustnessTemporalLogic(phi.lhs,R,time);

        cnt = 1; r = []; ind = []; alpha = {};

        while time(cnt) + phi.to < time(end)

            index = find(time >= phi.from & time <= phi.to);
            [rTmp,indexTmp] = max(r_(index));
            index = index(indexTmp);

            r = [r;rTmp]; ind = [ind;ind_(index)]; 
            alpha = [alpha,alpha_(index)];

            cnt = cnt + 1;
        end

    elseif strcmp(phi.type,'globally')

        [r_,alpha_,ind_] = robustnessTemporalLogic(phi.lhs,R,time);

        cnt = 1; r = []; ind = []; alpha = {};

        while time(cnt) + phi.to < time(end)

            index = find(time >= phi.from & time <= phi.to);
            [rTmp,indexTmp] = min(r_(index));
            index = index(indexTmp);

            r = [r;rTmp]; ind = [ind;ind_(index)]; 
            alpha = [alpha,alpha_(index)];

            cnt = cnt + 1;
        end

    elseif strcmp(phi.type,'until')

        [r1,alpha1,ind1] = robustnessTemporalLogic(phi.lhs,R,time);
        [r2,alpha2,ind2] = robustnessTemporalLogic(phi.rhs,R,time);

        
        cnt = 1; r = []; ind = []; alpha = {};
    
        while time(cnt) + phi.to < time(end)
    
            r = [r;-inf]; ind = [ind;0]; alpha = [alpha,{0}];
    
            index = find(time >= phi.from & time <= phi.to);
    
            for i = index'
                
                [r_,index_] = min(r1(1:i));
    
                if min(r_,r2(i)) > r(end)
                    if r_ < r2(i)
                        r(end) = r_; alpha{end} = alpha1{index_}; 
                        ind(end) = ind1(index_);
                    else
                        r(end) = r2(i); alpha{end} = alpha2{index_}; 
                        ind(end) = ind2(index_);
                    end
                end
            end
    
            cnt = cnt + 1;
        end
    end
end

function list = safe2unsafe(sets)
% convert a safe set defined by the union of multiple polytopes to an
% equivalent union of unsafe sets

    list = reverseHalfspaceConstraints(sets{1});

    for i = 2:length(sets)

        tmp = reverseHalfspaceConstraints(sets{i});

        list_ = {};

        for j = 1:length(tmp)
            for k = 1:length(list)
                if isIntersecting(list{k},tmp{j})
                    list_{end+1} = list{k} & tmp{j};
                end
            end
        end

        list = list_;
    end
end

function res = reverseHalfspaceConstraints(poly)
% get a list of reversed halfspace constraints for a given polytope

    res = {};
    poly = mptPolytope(poly);

    for i = 1:length(poly.P.b)
        res{end+1} = mptPolytope(-poly.P.A(i,:),-poly.P.b(i));
    end
end