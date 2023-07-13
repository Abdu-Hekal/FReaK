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
model.nResets=5;
%     model.solver.dt=10;
model.cp=[10,10];

model.inputInterpolation='previous';


% autokoopman settings
model.ak.obsType="rff";
model.ak.nObs=20;
model.ak.gridSlices=5;
model.ak.opt="grid"; %grid
model.ak.rank=[1,20,4];


%default optimizer options
solver = 'gurobi';  % gurobi, cplex, glpk
timeLimit = 60; %2000;
gapLimit = 1e-4; %0.1;
gapAbsLimit = 1e-10; %0.1;
solnLimit = Inf;
verb = 2;
model.solver.opts = sdpsettings('verbose', verb,'solver', solver, ...
    'gurobi.TimeLimit', timeLimit, ...
    'gurobi.MIPGap', gapLimit, ...
    'gurobi.MIPGapAbs', gapAbsLimit, ...
    'gurobi.SolutionLimit', solnLimit,...
    'gurobi.Method',3,...
    'gurobi.MIPFocus',3,...
    'gurobi.NumericFocus',3,...
    'usex0', 0 ...
    );

end