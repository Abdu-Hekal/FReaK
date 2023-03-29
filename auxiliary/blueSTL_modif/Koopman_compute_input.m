function [Sys, status] = Koopman_compute_input(Sys, controller)
% STLC_compute_input
%
% Input:
%       Sys: an STLC_lti instance
%       controller: a YALMIP optimizer object representing the system's
%                   optimization problem
%
% Output:
%       Sys: modified with additional model_data
%       params: controller data
%
% :copyright: TBD
% :license: TBD

%% Stuff
nu=Sys.nu;
nx=Sys.nx;

%% Time
ts=Sys.ts; % sampling time
L=Sys.L;  % horizon (# of steps)
t_model = floor(Sys.time(Sys.system_data.time_index)/Sys.ts)+1;
Sys.model_data.time = ((0:L)+max(t_model-L,0))*ts;

%% call solver
[sol_control, errorflag1] = controller{{}, []};
if(errorflag1==0)
    Sys.model_data.alpha = double(sol_control{1});
    Sys.model_data.X = double(sol_control{2});
    Sys.model_data.rob = double(sol_control{3});
elseif (errorflag1==1 || errorflag1==15||errorflag1==12)  % some error, infeasibility or else
    disp(['Yalmip error (disturbance too bad ?): ' yalmiperror(errorflag1)]); % probably there is no controller for this w
else
    disp(['Yalmip error: ' yalmiperror(errorflag1)]); % some other error
end
status = errorflag1;



end


