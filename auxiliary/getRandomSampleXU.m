
function [x0,u] = getRandomSampleXU(kfModel)
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
            cpVal = randPoint(kfModel.U(k),cp)';
            if all_steps > kfModel.cp(k)
                step = all_steps/kfModel.cp(k);
                assert(floor(step)==step,'number of control points (cp) must be a factor of T/ak.dt');
                u(:,k) = interp1((0:kfModel.ak.dt*step:kfModel.T-kfModel.ak.dt)', cpVal, linspace(0,kfModel.T-kfModel.ak.dt,all_steps)',kfModel.inputInterpolation,"extrap");
            else
                u(:,k) = cpVal;
            end
        end
    end
    u = [linspace(0,kfModel.T-kfModel.ak.dt,all_steps)',u];
else
    u = [];
end
end
