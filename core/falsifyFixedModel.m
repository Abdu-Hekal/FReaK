function [x0,u, kfModel] = falsifyFixedModel(A,B,g,kfModel)
% determine the most critical initial point and input trajectory of a
% Koopman linearized model for the given trajectory
%
% Input arguments:
%
%   -A:        system matrix for the Koopman linearized model
%   -B:        input matrix for the Koopman linearized model
%   -g:        function handle to the observable function for the Koopman
%              linearized model
%   -model:    koopman falsification model
%
% Output arguments:
%
%   -x0:       initial state for the most critical trajectory
%   -u:        piecewise constant inputs for the most critical trajectory

% compute reachable set for Koopman linearized model
R = reachKoopman(A,B,g,kfModel);
% determine most critical reachable set and specification
kfModel = mostCriticalReachSet(R,kfModel);

% extract most critical initial state and input signal
[x0,u] = falsifyingTrajectory(kfModel);

% %modification to test (delete me)
x = g(x0);
for i = 1:size(u,2)
    x = [x, A*x(:,end) + B*u(:,i)];
end
figure; hold on; box on;
for i=1:size(u,2)
    plot(R.zono{i})
end
plot(x(1,1:400),x(2,1:400),'r','LineWidth',2);
drawnow

end

% Auxiliary Functions -----------------------------------------------------
function R = reachKoopman(A,B,g,kfModel)
% compute reachable set for Koopman linearized model

%setup
dt=kfModel.dt;
R0=kfModel.R0;
U=kfModel.U;
tFinal=kfModel.T;
cpBool=kfModel.cpBool;

% compute initial set using Taylor model arithmetic
n = dim(R0);
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
R0 = polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:));

% compute reachable set
t = 0:dt:ceil(tFinal/dt)*dt;

set = cell(length(t),1); time = set; zono = set;

set{1} = R0; time{1} = interval(-dt/2,dt/2);

for i = 1:length(t)-1
    % AH edit to check if system has external input
    if B
        if kfModel.pulseInput
            cp_U = U.*cpBool(i,:)';
        else
            cp_U=U;
        end
        set{i+1} = A*set{i} + B*zonotope(cp_U);
    else
        set{i+1} = A*set{i};
    end
    time{i+1} = time{i} + dt;
end

% compute ouput set
for i = 1:length(set)
    set{i} = project(set{i},1:n);
    zono{i} = zonotope(set{i});

end

R.set = set; R.time = time; R.zono = zono;
end

function kfModel = mostCriticalReachSet(R,kfModel)
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
        maxStlSteps = min(maxStlTime(spec(i,1).set)/kfModel.dt,length(R.zono));
        %get prev solns
        prevSpecSol = kfModel.specSolns(kfModel.spec(i,1));
        %setup and run bluSTL
        tic
        %if no prev soln for this spec, setup alpha & stl milp vars/constrs
        try
            Sys=prevSpecSol.lti; %get previously setup milp problem with stl
            %setup offset and offset count
            if Sys.offsetCount==prevSpecSol.offsetCount %same subset to offset, then concatanate offsets
                Sys.offset=Sys.offset+prevSpecSol.realRob;
            else
            Sys.offset=prevSpecSol.realRob;
            end
            Sys.offsetCount=prevSpecSol.offsetCount;
        catch
            Sys=Koopman_lti(R.zono(1:maxStlSteps),kfModel.dt);
            if ~kfModel.pulseInput %if not pulse input, set cpBool
                Sys.cpBool=kfModel.cpBool;
            end
            Sys = setupAlpha(Sys);
            %convert disjunct stl from CORA format to blustl
            disjSet = disjunctiveNormalForm(spec(i,1).set); %is this necassary?
            bluStl = coraBlustlConvert(disjSet); %convert from cora syntax to blustl
            Sys.stlList = {bluStl};
            Sys=setupStl(Sys); %encode stl using milp
        end
        specSoln.lti=Sys; %store lti object with alpha problem info

        Sys.reachZonos=R.zono(1:maxStlSteps); %update reach zonos with new
        Sys = setupReach(Sys);
        
        kfModel.soln.milpSetupTime = kfModel.soln.milpSetupTime+toc;
        tic
        Sys = setupOptimizer(Sys,kfModel.solverOpts);
        Sys=optimize(Sys);
        kfModel.soln.milpSolvTime =kfModel.soln.milpSolvTime+toc;

        %get results
        rob_ = value(Sys.Pstl);
        alpha = value(Sys.Alpha);

        %clear solution from yalmip and assign only alpha for warmstarting
        %         yalmip('clearsolution')
        %         assign(Sys.Alpha,alpha)

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

function [x0,u] = falsifyingTrajectory(kfModel)
% extract the most critical initial state and input signal from the most
% critical reachable set and specification

%setup
R0=kfModel.R0;
U=kfModel.U;
cpBool=kfModel.cpBool;
set=kfModel.soln.set;
alpha=kfModel.soln.alpha;

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
%AH modification, checks if R0 is not exact, to avoid error
if isempty(generators(R0))
    x0 = center(R0);
else
    x0 = center(R0) + generators(R0)*alphaInit;
end

%check if kfModel has control input
if ~isempty(U)
    % determine most ctritical control input
    if ~isempty(set.Grest)

        if kfModel.pulseInput %if pulse input
            %initialise alpha to cpBool
            alphaU = reshape(cpBool,[],1);
            %input alpha returned by milp optimizernon
            allAlpha = alpha(size(set.G,2)+1:end);
            %find all nonzero elements in the relevant time horizon
            nonzero = find(alphaU,length(allAlpha));
            alphaU(nonzero) = allAlpha;
            %set all nonrelevant inputs after to zero
            alphaU(nonzero(end)+1:end)=0;
        else
            alphaU = alpha(size(set.G,2)+1:end);
        end

        U = zonotope(U); c_u = center(U); G_u = generators(U);
        alphaU = reshape(alphaU,[size(G_u,2),length(alphaU)/size(G_u,2)]);

        u = c_u + G_u*alphaU;
    else
        u = center(U);
    end
else
    u = [];
end
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
            maxTime=maxTime(obj.lhs)
        end
        if ~isempty(obj.rhs)
            maxTime=maxTime(obj.rhs)
        end
    end

end
