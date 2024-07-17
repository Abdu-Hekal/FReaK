function [tout, yout, uout, simTime] = correctiveSim(obj, koopModel, xTarget)
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
%    KoopModel  - learnt koopman model with A and B matrices and observables.
%    xTarget   - target trajectory to follow.
%
% Outputs:
%    tout - Time vector of simulation
%    yout - Output vector of simulation
%    uout - input vector of simulation (after corrective feedback)
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

simTime=0;%intialize simulation time
tout=[]; %empty list to store time points
yout=[]; %empty list to store trajectory
uout=[]; %empty list to store inputs
n = obj.R0.dim; %number of variables
origT=obj.T; %store original time horizon for problem

t0=0; %start time
error=zeros(1,n); %initial error between real traj and target
x0=xTarget(:,1); %initial x
xFinal=x0(1:n);
%loop over all time points in trajectory
for ii=2:size(xTarget,2)
    x1 = xTarget(:,ii)+koopModel.g((error)');

    x1x0 =(x1 - koopModel.A* x0);
    u = pinv(koopModel.B(1:n,:)) * x1x0(1:n); 
    u = max(obj.U.inf,min(obj.U.sup,u)); %ensure that extrapolation is within input bounds

    if ii==2 %first iteration of loop
        obj.T=t0+obj.ak.dt-obj.dt; %modify time horizon for small steps (minus last step)
        uTimes=(t0:obj.dt:obj.T)';
        u=[uTimes,repmat(u',numel(uTimes),1)];
        uout=u(1,:); %autokoopman input list
    elseif ii==size(xTarget,2) %last iteration of loop
        obj.T=t0+obj.ak.dt+obj.dt; %modify time horizon for small steps
        uTimes=(t0:obj.dt:obj.T)';
        %add previous inp as first step, then add new input
        u=[uTimes, [uout(end,2:end);repmat(u',numel(uTimes)-1,1)]];
        uout=[uout;u(2,:);u(end,:)]; %autokoopman input list
    else
        obj.T=t0+obj.ak.dt; %modify time horizon for small steps
        uTimes=(t0:obj.dt:obj.T)';
        %add previous inp as first step, then add new input
        u=[uTimes, [uout(end,2:end);repmat(u',numel(uTimes)-1,1)]];
        uout=[uout;u(2,:)]; %autokoopman input list
    end
    [tout_, yout_, simTime_,xFinal]=simulate(obj, xFinal, u, t0);

    error= x1(1:n)'-yout_(end,:); %difference between target point and real simulated point

    %note: we remove last element to avoid duplication with new simulation
    tout=[tout(1:end-1);tout_]; %update list of all time points simulated
    yout=[yout(1:end-1,:);yout_]; %update list of real trajectory
%     u=[t0,u(1,2:end)];
%     uout=[uout(1:end-1,:);u];
    
    t0=tout(end); %new time point is last time point in current simulation
    x0=yout_(end,:)'; %current point is updated to new point
    x0=koopModel.g(x0);
    simTime=simTime+simTime_; %update simulation time
end

obj.T=origT; %revert KF property 'T' back to full time horizon for problem

end