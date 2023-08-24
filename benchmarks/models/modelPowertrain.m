function model = modelPowertrain()
% modelCars - model parameters for the Fuel Control of an Automotive Powertrain (AFC).
%
% Syntax:
%       model = modelNeural()
%
% Description:
%       Model parameters for the AFC benchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model
%
%------------------------------------------------------------------

model = KF_model(@run_powertrain);
model.R0 = interval(0,0);
model.U = interval([0;900],[61.1;1100]);

model.T=50;
model.dt = 0.01;
model.ak.dt= 5;
%     model.solver.dt=10;
model.cp=[10,10];

model.inputInterpolation='previous';

end