function controller = reach_get_controller(Sys,~)
% STLC_get_controller constructs the controller object for an STLC_lti instance
%
% Input:
%       Sys: an STLC_lti instance
%
% Output:
%       controller: a YALMIP optimizer object that solves the STL-constrained
%                   optimal control problem for Sys
%
% :copyright: TBD
% :license: TBD


%% Time
L=Sys.L;   % horizon (# of steps)
ts=Sys.ts; % sampling time

%% System dimensions and variables
nx=size(Sys.x0,1);
% variables
X = sdpvar(nx, L+1);
Alpha = sdpvar(1, size(Sys.reach_zonos{end}.generators,2));

%% STL formula
Fstl = [];
Pstl = [];
var = struct('X',X);

stl_list= STLC_parse_stl_labels(Sys);
M = Sys.bigM;


for i = 1:numel(stl_list)
    phi = STLformula('phi', stl_list{i});

    %what time steps to compute robustness for? does this work? why?
    [Fphi, Pphi] = Koopman_MILP_robust(phi, 1, L+1, ts, var,M);
    
    Pstl = [Pstl; Pphi];
    Fstl = [Fstl Fphi];

end


%% Reachset constraints
Freach = [];

%constraints for Alpha
for k=1:size(Alpha,2)
    Freach = [Freach, -1 <= Alpha(k) <= 1];
end

% Constraints for reachable set
for k=1:L+1
    % x = c + G * \alpha, 
    c = Sys.reach_zonos{k}.center;
    G = Sys.reach_zonos{k}.generators;

    Freach = [Freach, X(:,k) == (c+G*Alpha(1:size(G,2))')];
end

options = Sys.solver_options;
output_controller =  {Alpha,X,Pstl};


if numel(stl_list) == 0
    Pstl = sdpvar(1,1);
end

%% Objective function
obj = sum(sum(Pstl(:,1:end))); 

controller = optimizer([Freach, Fstl],obj,options,[], output_controller);


