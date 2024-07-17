function model = modelSC2()
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
    
    model = KF(@runSteamcondenser);
    model.R0 = interval([80;107.9;9062.6;90],[80;107.9;9062.6;90]); 
    model.U = interval(3.99,4.01); 

    model.T=35; 
    model.ak.dt= 0.1; %0.25;
%     model.solver.dt=10;
    model.cp=20;

end