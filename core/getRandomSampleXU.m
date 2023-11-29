function [x0,u] = getRandomSampleXU(kfModel)
% getRandomSampleXU - Generate random initial state and input for the KF model
%
% Syntax:
%    [x0, u] = getRandomSampleXU(kfModel)
%
% Description:
%    This function generates a random initial state (x0) and an associated
%    input signal (u) for the Koopman model based on the provided KF object
%    and its parameters.
%
% Inputs:
%    kfModel - KF object containing the Koopman model and various
%              parameters needed for the falsification process.
%
% Outputs:
%    x0 - Randomly generated initial state.
%    u - Randomly generated input signal associated with the initial state.
%
% Author:      Abdelrahman Hekal
% Written:     28-February-2023
% Last update: [Date]
% Last revision: [Date]
%
% ------------- BEGIN CODE --------------

%generate random initial set
x0 = randPoint(kfModel.R0);
%generate random input if kfModel has input.
u=[];
if ~isempty(kfModel.U)
    all_steps = kfModel.T/kfModel.ak.dt;
    if kfModel.pulseInput
        u = randPoint(kfModel.U,all_steps)';
        u = u.*kfModel.cpBool;
    else %piecewise constant input
        for k=1:length(kfModel.cp)
            cp = min(all_steps, kfModel.cp(k)); %control points is minimum of maximum control points and koopman time points (can't have more control points than steps)
            cpVal = randPoint(kfModel.U(k),cp+1)'; %add 1 to cp for last timestep.
            if all_steps > kfModel.cp(k)
                step = all_steps/kfModel.cp(k);
                assert(floor(step)==step,'number of control points (cp) must be a factor of T/ak.dt');
                u(:,k) = interp1((0:kfModel.ak.dt*step:kfModel.T)', cpVal, linspace(0,kfModel.T,all_steps+1)',kfModel.inputInterpolation,"extrap");
            else
                u(:,k) = cpVal;
            end
        end
    end
    u = [linspace(0,kfModel.T,all_steps+1)',u];
else
    u = [];
end
end
