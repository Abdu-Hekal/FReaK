function model = modelCars()
% modelCars - model parameters for the chasing cars benchmark
%
% Syntax:
%       model = modelCars()
%
% Description:
%       Model parameters for the chasing cars benchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model      
%
%------------------------------------------------------------------
    
    model = KF_model(@run_cars); 
    model.R0 = interval([0;10;20;30;40],[0;10;20;30;40]); 
    model.U = interval([0;0],[1;1]); 

    model.T=100; 
    model.dt = 0.01; 
    model.ak.dt=2.5; %5
%     model.solver.dt=10;
    model.cp=[10000 10000];

    x = stl('x',5); 
    eq = globally(x(5)-x(4)<=40,interval(0,100));
    model.spec = specification(eq,'logic');

end