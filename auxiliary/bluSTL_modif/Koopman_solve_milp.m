function [model_data, status] = Koopman_solve_milp(Sys,milp)
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

%setup parameters
reach_zonos=Sys.reach_zonos;
c = zeros(length(reach_zonos{end}.center),length(reach_zonos));
for k=1:length(reach_zonos)
    c(:,k)=reach_zonos{k}.center;
end
G=reach_zonos{end}.generators;

%% call solver
[sol_control, errorflag1] = milp{{c,G}};
if(errorflag1==0)
    model_data.alpha = double(sol_control{1});
    model_data.X = double(sol_control{2});
    model_data.rob = double(sol_control{3});
else
    disp(['Yalmip error: ' yalmiperror(errorflag1)]); % some other error
end
status = errorflag1;

end


