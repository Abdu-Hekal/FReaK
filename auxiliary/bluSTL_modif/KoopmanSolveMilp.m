function [model_data, status] = KoopmanSolveMilp(milp)
% KoopmanSolveMilp
%
% Input:
%       Sys: a Koopman_lti instance
%       milp: a YALMIP optimizer object representing the system's
%                   optimization problem
%
% Output:
%       model_data: struct with results
%       status: yalmip error flag
%
% :copyright: TBD
% :license: TBD

%% call solver
[sol_control, errorflag1] = milp{{}, []};
if(errorflag1==0)
    model_data.alpha = double(sol_control{1});
    model_data.X = double(sol_control{2});
    model_data.rob = double(sol_control{3});
else
    disp(['Yalmip error: ' yalmiperror(errorflag1)]); % some other error
end
status = errorflag1;
end


