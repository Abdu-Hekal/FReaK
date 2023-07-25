function model = modelSC()
% modelCars - model parameters for the Steam condenser with Recurrent Neural Network Controller (SC).
%
% Syntax:
%       model = modelSC()
%
% Description:
%       Model parameters for the Steam condenserbenchmark.
%
% Output Arguments:
%
%       -model:             a koopman falsification model      
%
%------------------------------------------------------------------
    
    model = KF_model(@run_steamcondenser);
    model.R0 = interval([80;107.9;9062.6;90],[80;107.9;9062.6;90]); 
    model.U = interval(3.99,4.01); 

    model.T=35; 
    model.dt = 0.01; 
    model.ak.dt= 2.5; %0.25;
    model.nResets=5;
%     model.solver.dt=10;
    model.cp=3500;

end