function model = modelVanderpol()
% model_vanderpol - model parameters for the vanderpol benchmark
%
% Syntax:
%       model = model_vanderpol()
%
% Description:
%       Model parameters for the vanderpol benchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model
%
%------------------------------------------------------------------

    model = KF('vanderpol');
    model.R0 = interval([1.25;2.25],[1.55;2.35]); 
    
    model.T=7; 
    model.dt = 0.01;
    model.ak.dt=0.1;

    model.spec = specification(halfspace([-1;0],-2.095),'unsafeSet');
    x = stl('x',2);
    eq = globally(x(1) < 2.095,interval(0,7));
    model.spec = specification(eq,'logic');

end