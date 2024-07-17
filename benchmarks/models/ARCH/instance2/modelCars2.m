function model = modelCars2()
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
    
    model = KF('cars'); 
    model.R0 = interval([0;10;20;30;40],[0;10;20;30;40]); 
    model.U = interval([0;0],[1;1]); 

    model.T=100; 
    model.ak.dt=1; %5
    model.cp=[20 20];

    x = stl('x',5); 
    eq = globally(x(5)-x(4)<=40,interval(0,100));
    model.spec = specification(eq,'logic');

end