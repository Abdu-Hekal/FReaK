function model = modelPacemaker()
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
    
    model = KF(@runPacemaker); 
    model.R0 = interval(0,0); 
    model.U = interval(50,90); 

    model.T=10; 
    model.ak.dt=0.1;
    model.cp=1000;
    model.inputInterpolation='previous';

    x = stl('x',1); 
    eq = globally(x(1)<=15,interval(0,10)) & finally(x(1)>=8,interval(0,10));
    model.spec = specification(eq,'logic');

%     model.solver.autoAddTimePoints=true;
%     model.solver.autoAddConstraints=2;

    model.verb=3;

end