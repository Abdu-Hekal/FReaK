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
%       -model:             a koopman falsification model containing
%                           following properties
%
%           -.name:         name of simulink model
%           -.T:            time horizon for simulation
%           -.R0:           initial set (CORA class interval)
%           -.U:            set of admissible control inputs (class:
%                           interval or Zonotope)
%           -.N:            no. of (piecewise constant) inputs, 
%                           set to 1 if no inputs
%           -spec:          specification defined as an object of the CORA specification
%                           class (includes safe sets, unsafe sets, and temporal logic)
%
%------------------------------------------------------------------
    
    % initalize training data and model specific information
    % if U don't exist for the model, set as "interval(0,0)"
    % N: no. of (piecewise constant) inputs, set to 1 if no inputs
    % spec: specification of (unsafe) set
    model.name = 'vanderpol'; %name of the simulink model
    model.T=7; 
    model.R0 = interval([1.25;2.25],[1.55;2.35]); 
    model.U = interval(0,0); 
    model.N=1;
    model.spec = specification(halfspace([-1;0],-2.095),'unsafeSet');

end