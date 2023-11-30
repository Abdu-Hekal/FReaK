function [obj, koopModel] = learnKoopModel(obj, trainset)
% learnKoopModel - Symbolically compute the Koopman linearized model and
%   observables function using Autokoopman (python library)
%
% Syntax:
%    [obj, koopModel] = learnKoopModel(obj, trainset)
%
% Description:
%    This function symbolically computes the Koopman linearized model and
%    observables using the Random Fourier Features (RFF) method. It sets up
%    the AutoKoopman settings, runs the AutoKoopman Python script, and loads
%    the computed A, B, u, and w matrices. It then creates and saves the
%    observables function, returning the updated obj, and a struct koopModel
%    with A, B, and the observables function g.
%
% Inputs:
%    obj - KF object containing Koopman model and various parameters
%              needed for the symbolic RFF process.
%    trainset - Training set used to compute the Koopman model.
%
% Outputs:
%    obj - Updated KF object.
%    koopModel - struct representing koopman model with following:
%       A - State transition matrix of the Koopman linearized model.
%       B - Input matrix of the Koopman linearized model.
%       g - observables function of the Koopman linearized model.
%
% Example:
%    [obj, A, B, g] = learnKoopModel(obj, trainset);
%
% See also: falsify
%
% Author: Abdelrahman Hekal
% Written: 28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % AutoKoopman settings and run
    if ~isempty(obj.U) %check if obj has inputs
        inputs_list = trainset.XU;
    else
        inputs_list = string(missing);
    end 
    tic
    pyrunfile("run_autokoopman.py",'koopman_model',times=trainset.t,trajectories=trainset.X, param_dict=obj.ak,inputs_list=inputs_list);
    obj.soln.koopTime= obj.soln.koopTime+toc;
    %TODO: remove need for creating file by returning A,B,u,w functions
    %from python file.
    load("autokoopman_model.mat", "A","B","u","w")

    % create observables function
    n = size(obj.R0,1); %number of variables
    g_ = @(x) sqrt(2/obj.ak.nObs)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    g = matlabFunction(g,'Vars',{xSym});

    %store koopman model
    koopModel.A=A; koopModel.B=B; koopModel.g=g;
end