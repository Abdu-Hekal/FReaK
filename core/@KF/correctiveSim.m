function [tout, yout, simTime] = correctiveSim(obj, koopModel, xTarget)
% correctiveSim - Simulate the model associated with a Koopman Falsification
% object using corrective feedback.
%
% Syntax:
%    [tout, yout, simTime] = correctiveSim(obj, koopModel, xTarget)
%
% Description:
%    This function simulates the model associated with a Koopman
%    Falsification (KF) object using corrective feedback. Given a KF
%    object, a koopman model and a target trajectory, the simulator tries
%    to follow the trajectory as best as possible by employing a corrective
%    feedback loop which adjusts inputs after each step   
%    The simulated model can be either a Simulink model
%    or a custom function handle. Note that a custom function can be
%    used for simulation by passing the function handle. Ensure that the
%    outputs are consistent with the simulation method used for the
%    specific model.
%
% Inputs:
%    KoopModel  - learnt koopman model with A and B matrices.
%    xTarget   - target trajectory to follow.
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    simTime  - time taken for simulation
%
%
% See also: falsify, simulate
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: 4-December-2023
% Last revision: ---
%------------- BEGIN CODE --------------

tsim = (0:obj.dt:obj.T)'; %time points for interpolating trajectory
xTarget = interp1(xTarget(:,1),xTarget(:,2:end),tsim,"linear"); %interpolate trajectory

for 

end