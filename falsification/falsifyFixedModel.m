function [x0,u] = falsifyFixedModel(A,B,g,spec,R0,U,tFinal,dt)

    robustnessPolytope(1,2);

    % compute reachable set for Koopman linearized model
    R = reachKoopman(A,B,g,R0,U,tFinal,dt);

    % determine most critical reachable set and specification
    [set,spec] = mostCriticalReachSet(R,spec);

    % extract most critical initial state and input signal
    [x0,u] = falsifyingTrajectory(R0,U,set,spec);

    % debug
    figure; hold on; box on;
    for i = 1:length(R.zono)
        plot(R.zono{i},[1,2],'b');
    end
    plot(set,[1,2],'g');
    x = [g(x0)];
    for i = 1:length(R.zono)
        x = [x,A*x(:,end)];
    end
    plot(x(1,:),x(2,:),'r');

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

function [setCrit,specCrit] = mostCriticalReachSet(R,spec)
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
                    dist = infimum(interval(set.P.A*R.zono{j})) - set.P.b;
                    [rob_,ind] = max(dist);

                    % check if smaller than current robustness
                    if rob_ < rob
                        specCrit = halfspace(set.P.A(ind,:),set.P.b(ind));
                        setCrit = R.set{j};
                        rob = rob_;
                    end
                end
            end

        elseif strcmp(spec(i,1).type,'safeSet')


        elseif strcmp(spec(i,1).type,'logic')

        else
            error('This type of specification is not supported!');
        end
    end
end

function [x0,u] = falsifyingTrajectory(R0,U,set,spec)
% extract the most critical initial state and input signal from the most
% critical reachable set and specification

    % determine most critical initial state
    temp = spec.c'*set.G;
    alpha = zeros(size(set.expMat,1),1);

    ind = find(sum(set.expMat,1) == 1);

    for i = 1:size(set.expMat,1)
        for j = ind
            if set.expMat(i,j) == 1
                alpha(i) = -sign(temp(j));
            end
        end
    end

    R0 = zonotope(R0);
    x0 = center(R0) + generators(R0)*alpha;

    % determine most ctritical control input
    if ~isempty(set.Grest)
        alpha = -sign(spec.c' * set.Grest);
    
        U = zonotope(U); c_u = center(U); G_u = generators(U);
        alpha = reshape(alpha,[size(G_u,2),length(alpha)/size(G_u,2)]);

        u = c_u + G_u*alpha;
    else
        u = center(U);
    end
end

function [r,alpha] = robustnessPolytope(P,Z)
% compute robustness of the zonotope Z with respect to an unsafe polytope P
    
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