function [koopModel,koopTime] = learnKoopModel(obj, trainset)
% learnKoopModel - Symbolically compute the Koopman linearized model and
%   observables function using Autokoopman (python library)
%
% Syntax:
%    [koopModel,koopTime] = learnKoopModel(obj, trainset)
%
% Description:
%    This function symbolically computes the Koopman linearized model and
%    observables using the Random Fourier Features (RFF) method. It sets up
%    the AutoKoopman settings, runs the AutoKoopman Python script, and loads
%    the computed A, B, u, and w matrices. It then creates and saves the
%    observables function, a struct koopModel with A, B, and the observables
%    function g. It also returns the time taken for autokoopman to run.
%
% Inputs:
%    obj - KF object containing Koopman model and various parameters
%              needed for the symbolic RFF process.
%    trainset - Training set used to compute the Koopman model.
%
% Outputs:
%    koopModel - struct representing koopman model with following:
%       A - State transition matrix of the Koopman linearized model.
%       B - Input matrix of the Koopman linearized model.
%       g - observables function of the Koopman linearized model.
%    koopTime - time taken for autokoopman to tune and learn model
%
% Example:
%    [koopModel,koopTime] = learnKoopModel(obj, trainset)
%
% See also: falsify
%
% Author: Abdelrahman Hekal
% Written: 28-February-2023
% Last update: 4-December-2023
% Last revision: ---

%------------- BEGIN CODE --------------

% AutoKoopman settings and run
if ~isempty(obj.U) %check if obj has inputs
    inputs_list = trainset.XU;
else
    inputs_list = string(missing);
end
if obj.ak.weighted
    robList = trainset.Rob;
    stateWeights=ones(dim(obj.R0),1)*numel(trainset.t{1});
    % if only one spec and it is stl and state weights options set to true, then use weight states
    % according to variable hierarchy
    if numel(obj.spec)==1 && strcmp(obj.spec(1).type,'logic') && obj.ak.stateWeighted
        stateWeights=getVarHierarchy(obj.spec(1).set,obj.ak.dt,obj.T,1,stateWeights);
    end
else
    robList = string(missing);
    stateWeights = string(missing);
end
tic
pyrunfile("run_autokoopman.py",'koopman_model',times=trainset.t,trajectories=trainset.X,param_dict=obj.ak,inputs_list=inputs_list,rob_list=robList,state_weights=stateWeights);
koopTime=toc;
%TODO: remove need for creating file by returning A,B,u,w functions
%from python file.
load("autokoopman_model.mat", "A","B","u","w")

% create observables function
n = dim(obj.R0); %number of variables
g_ = @(x) sqrt(2/obj.ak.nObs)*cos(w*x + u');
g = @(x) [x; g_(x)];
xSym = sym('x',[n,1]);
g = g(xSym);

% save observables in function
g = matlabFunction(g,'Vars',{xSym});

%store koopman model
koopModel.A=A; koopModel.B=B; koopModel.g=g;
end

function stateWeights=getVarHierarchy(set,dt,T,count,stateWeights)
if isempty(set) || isnumeric(set)
    return;
elseif strcmp(set.type,'variable')
    %assign state weights for variables based on predecessor nodes
    idx=strcmp(set.var,getVariables(set));
    stateWeights(idx)=count;
else
    %temporal property
    if ~isempty(set.from) && ~isempty(set.to)
        %ensure time limits of temporal operators are within time horizon [0,T]
        if set.from<0
            set.from=0;
        end
        if set.to>T
            set.to=T;
        end
        %increase count based on number of leaves
        count = count * ceil((set.to-set.from)/dt);
    end
    stateWeights=getVarHierarchy(set.lhs,dt,T,count,stateWeights);
    stateWeights=getVarHierarchy(set.rhs,dt,T,count,stateWeights);
end
end