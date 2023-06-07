function model = model_vanderpol()
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

    model = KF_model('vanderpol');
    model.R0 = interval([1.25;2.25],[1.55;2.35]); 
    
    model.T=7; 
    model.dt = 0.01;

    model.spec = specification(halfspace([-1;0],-2.095),'unsafeSet');

end