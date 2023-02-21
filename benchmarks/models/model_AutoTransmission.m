function model = model_AutoTransmission()
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
%       -model:             a structure containing following options
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
    model.name = 'Autotrans_shift'; %name of the simulink model
    model.T=30; 
    model.R0 = interval([0;1000;1],[0;1000;1]); 
    model.U = interval([0;0],[100;320]); 
    model.N=3000;

    x = stl('x',3);
    eq = globally(x(2) < 4750,interval(0,30));
    model.spec = specification(eq,'logic');
%     model.spec = specification(halfspace([0 -1 0],-4750),'unsafeSet');

end