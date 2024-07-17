function [x0,u] = falsifyingTrajectory(obj,soln)
% falsifyingTrajectory - Extract the most critical initial state and input
%   signal from the set of most critical alpha values of a reachable set
%
% Syntax:
%    [x0, u] = falsifyingTrajectory(obj,soln)
%
% Description:
%    This function is responsible for extracting the most critical initial state
%    (x0) and the associated control input signal (u) from the most critical
%    reachable set and specification obtained during the falsification process.
%    It is a key step in the core falsification procedure.
%
% Inputs:
%    obj - KF object containing the Koopman model and various 
%       parameters needed for the falsification process.
%    soln - soln struct with critical set, alpha and/or input
%
% Outputs:
%    x0 - Most critical initial state.
%    u - Control input signal associated with the most critical trajectory.
%
% Example:
%    [x0, u] = falsifyingTrajectory(obj,soln);
%
% See also: falsify
%
% Author:      Niklas Kochdumper, Abdelrahman Hekal
% Written:     28-February-2023
% Last update:  4-Decemeber-2023
% Last revision:---

%------------- BEGIN CODE --------------

%setup
R0=obj.R0;
U=obj.U;
cpBool=obj.cpBool;
set=soln.set;
alpha=soln.alpha;

R0 = zonotope(R0);
%checks if R0 is exact
if isempty(generators(R0))
    x0 = center(R0);
else
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
    G_R0_ = generators(R0);
    %FIXME: AH modification, ensures that generators is n*n matrix, appends zeros for
    %dimensions with exact x0
    G_R0 = zeros(size(generators(R0),1));
    [row ,col]=find(G_R0_);
    for ii=1:numel(row)
        G_R0(row(ii),row(ii)) = G_R0_(row(ii),col(ii));
    end
    x0 = center(R0) + G_R0*alphaInit;
end
%check if obj has control input
if ~isempty(U)
    if isempty(soln.u)
        % determine most ctritical control input
        if ~isempty(set.Grest)

            if obj.pulseInput %if pulse input
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
        % control input already provided by solver
        u=soln.u;
    end
else
    u = [];
end
if ~isempty(u)
    all_steps = obj.T/obj.ak.dt;
    %check that correct number of inputs is returned, sometimes we have a
    %final dummy input so that number of inputs is equal to state variables
    %for MILP stl encoding, we remove it.
    if size(u,2) == all_steps
    elseif size(u,2)==all_steps+1 %inputs returned include last time step
        u=u(:,1:end-1);
    else
        error('incorrect number of inputs returned from solver, investigate error')
    end
    tp_=linspace(0,obj.T-obj.ak.dt,all_steps); %time points without last time step
    tp = linspace(0,obj.T,all_steps+1);
    u = interp1(tp_',u',tp',obj.inputInterpolation,"extrap"); %interpolate and extrapolate input points
    u =  max(obj.U.inf',min(obj.U.sup',u)); %ensure that extrapolation is within input bounds
    u = [tp',u];
end
end