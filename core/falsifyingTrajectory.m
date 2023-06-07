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
%append zeros for remaining steps if no input given.
all_steps = kfModel.T/kfModel.ak.dt;
u = [u';zeros(size(u,1),all_steps-size(u,2))']; 
%append time points as first column
u = [linspace(0,kfModel.T-kfModel.ak.dt,all_steps)',u];
end