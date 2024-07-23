function R = reachKoopman(obj,koopModel,prevSoln)
% reachKoopman - Compute the reachable set for the Koopman linearized model
%
% Syntax:
%    R = reachKoopman(obj,koopModel)
%
% Description:
%    This function calculates the reachable set for the Koopman linearized
%    model, given the system matrices A and B, the observables function g,
%    and the Koopman Falsification (KF) model parameters. The reachable set
%    is computed using polynomial zonotopes and updated over discrete
%    time steps.
%
% Inputs:
%    obj - KF object containing various parameters needed for the
%              falsification process.
%    koopModel - struct representing koopman model with following:
%       A - State transition matrix of the Koopman linearized model.
%       B - Input matrix of the Koopman linearized model.
%       g - observables function of the Koopman linearized model.
%    prevSoln - A struct representing critical (best) soln from previous iteration
%    with following:
%       koopSolver: KoopSolver object used
%       rob: best robustness value found on the (approximate) Koopman model
%       alpha: alpha values for critical trajectory
%       u: critical input
%       x: critical Koopman trajectory
%       set: critical polynomial zonotope set
%
% Outputs:
%    R - Struct containing the reachable set information, including
%       polynomial zonotopes, time points, and zonotopes.
%
% Example:
%    R = reachKoopman(obj, koopModel);
%
% See also: falsify, critAlpha
%
% Author: Niklas Kochdumper, Abdelrahman Hekal
% Written: 28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%setup
A = koopModel.A;
B = koopModel.B;
g = koopModel.g;
dt=obj.ak.dt;
R0=obj.R0;
U=obj.U;
tFinal=obj.T;
cpBool=obj.cpBool;
tayOrder = obj.reach.tayOrder;

% compute initial set using Taylor model arithmetic
n = dim(R0); dig = length(num2str(n));
names = {}; for i = 1:n; names{i,1} = ['x',num2str(i,['%0',num2str(dig), '.f'])]; end

if all(rad(R0) == 0)
    c = g(center(R0));
    R0 = polyZonotope(c,zeros(length(c),1),[],zeros(n,1));
else
    tay = taylm(R0,tayOrder,names);
    tay = g(tay);
    R0 = polyZonotope(tay);
    R0 = polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:));
end

% check if B is empty anf U isn't, then set B to zeros to avoid errors in
% optimization
if isempty(B) && ~isempty(U)
    B=zeros(size(A,1),dim(U));
end

% compute reachable set
t = 0:dt:ceil(tFinal/dt)*dt;

set = cell(length(t),1); time = set; zono = set;

set{1} = R0; time{1} = interval(-dt/2,dt/2);

for i = 1:length(t)-1
    % check if system has external input
    if ~isempty(B)
        if obj.pulseInput
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
% get alpha factors for initial set
if isempty(prevSoln) || isempty(generators(zonotope(obj.R0))) 
    alphaInit=[];
else
    prevSet = prevSoln.set;
    alphaInit = zeros(size(prevSet.expMat,1),1);
    temp = prod(ones(size(prevSet.expMat))-mod(prevSet.expMat,2),1);
    expMat = prevSet.expMat(:,temp == 0);

    ind = find(sum(expMat,1) == 1);

    for i = 1:size(prevSet.expMat,1)
        for j = ind
            if expMat(i,j) == 1
                alphaInit(i) = prevSoln.alpha(j);
            end
        end
    end
end
% compute ouput set
for i = 1:length(set)
    set{i} = project(set{i},1:n);
    if ~isempty(alphaInit)
        zono{i} = zonotopeEnclosureMinLinErr(set{i},alphaInit,false);
        %         if i==1
        %             zono{i} = zonotopeEnclosureNew(set{i},alphaInit');
        %         else
        %             zono{i} = zonotope(set{i});
        %         end
    else
        zono{i} = zonotope(set{i});
    end
end

R.set = set; R.time = time; R.zono = zono;
end
