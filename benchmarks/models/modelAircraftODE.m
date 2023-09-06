function model = modelAircraftODE()
% modelModulator - model parameters for the aircraft ode benchmark
%
% Syntax:
%       model = modelModulator()
%
% Description:
%       Model parameters for the aircraft ode benchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model
%
%------------------------------------------------------------------

model = KF_model(@run_aircraft);
model.R0 = interval([200;-10;120],[260;10;150]);
model.U = interval([34386;0],[53973;16]);

model.T=4;
model.dt = 0.01;
model.ak.dt=0.1;
model.cp=[400 400];

x = stl('x',3);
eq = ~(globally(x(1) >=240 & x(1) <=250,interval(0,4)) & finally(x(1)>=240 & x(1)<=240.1,interval(3.5,4)));
%     eq = conjunctiveNormalForm(~(globally(x(1) >=240 & x(1) <=250,interval(0,4)) & finally(x(1)>=240 & x(1)<=240.1,interval(3.5,4))));

model.spec = specification(eq,'logic');


end