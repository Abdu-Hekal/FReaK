function model = modelModulator()
% modelModulator - model parameters for the modulator benchmark
%
% Syntax:
%       model = modelModulator()
%
% Description:
%       Model parameters for the modulator benchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model      
%
%------------------------------------------------------------------
    
    model = KF(@runModulator);
    model.R0 = interval([-0.1;-0.1;-0.1],[0.1;0.1;0.1]); 
    model.U = interval(-0.45,0.45); 

    model.T=9; 
    model.dt = 0.01;
    model.ak.dt=0.09; %2.5
    model.cp=100;

    x = stl('x',3);
    eq = globally(x(1) >=-1 & x(1) <=1 & x(2) >=-1 & x(2) <=1 & x(3) >=-1 & x(3) <=1,interval(0,9));

    model.spec = specification(eq,'logic');


end