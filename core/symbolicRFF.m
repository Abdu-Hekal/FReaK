function [kfModel, A, B,g] = symbolicRFF(kfModel, trainset)
% symbolicRFF - Symbolically compute the Koopman linearized model and
%   observables function using Autokoopman (python library)
%
% Syntax:
%    [kfModel, A, B, g] = symbolicRFF(kfModel, trainset)
%
% Description:
%    This function symbolically computes the Koopman linearized model and
%    observables using the Random Fourier Features (RFF) method. It sets up
%    the AutoKoopman settings, runs the AutoKoopman Python script, and loads
%    the computed A, B, u, and w matrices. It then creates and saves the
%    observables function, returning the updated kfModel, A, B, and the
%    observables function g.
%
% Inputs:
%    kfModel - KF object containing Koopman model and various parameters
%              needed for the symbolic RFF process.
%    trainset - Training set used to compute the Koopman model.
%
% Outputs:
%    kfModel - Updated KF object.
%    A - Koopman matrix A.
%    B - Input matrix B.
%    g - Observables function.
%
% Example:
%    [kfModel, A, B, g] = symbolicRFF(kfModel, trainset);
%
% See also: falsify
%
% Author: Abdelrahman Hekal
% Written: 28-February-2023
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % AutoKoopman settings and run
    if ~isempty(kfModel.U) %check if kfModel has inputs
        inputs_list = trainset.XU;
    else
        inputs_list = string(missing);
    end 
    tic
    pyrunfile("run_autokoopman.py",'koopman_model',times=trainset.t,trajectories=trainset.X, param_dict=kfModel.ak,inputs_list=inputs_list);
    kfModel.soln.koopTime= kfModel.soln.koopTime+toc;
    %TODO: remove need for creating file by returning A,B,u,w functions
    %from python file.
    load("autokoopman_model.mat", "A","B","u","w")

    % create observables function
    n = size(kfModel.R0,1); %number of variables
    g_ = @(x) sqrt(2/kfModel.ak.nObs)*cos(w*x + u');
    g = @(x) [x; g_(x)];
    xSym = sym('x',[n,1]);
    g = g(xSym);
    
    % save observables in function
    g = matlabFunction(g,'Vars',{xSym});
end