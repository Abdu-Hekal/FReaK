function milp = reach_setup_milp(Sys)
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


%% setup
L=Sys.L;   % horizon (# of steps)
dt=Sys.dt; % sampling time
cp_bool=Sys.cp_bool;
reach_zonos=Sys.reach_zonos;

%% System dimensions and variables
nx=Sys.nx;
% variables
X = sdpvar(nx, L+1);
Alpha = sdpvar(1, size(reach_zonos{end}.generators,2));

%% STL formula
Fstl = [];
Pstl = [];
var = struct('X',X);

stl_list= Koopman_parse_stl_labels(Sys);
M = Sys.bigM;


for i = 1:numel(stl_list)
    phi = STLformula('phi', stl_list{i});

    [Fphi, Pphi] = Koopman_MILP_robust(phi, 1, L+1, dt, var,M);
    
    Pstl = [Pstl; Pphi];
    Fstl = [Fstl Fphi];

end


%% Reachset constraints
Falpha = [];

%constraints for Alpha
for k=1:size(Alpha,2)
    Falpha = [Falpha, -1 <= Alpha(k) <= 1];
end

%constraint for control points
if ~isempty(cp_bool) %piecewise constant signal
    %first alpha corresponding to an input
    k = size(reach_zonos{1}.generators,2)+1;
    for col=1:size(cp_bool,2)
        for row=1:size(cp_bool,1)
            %if bool is zero constrain alpha to be same value as prev
            if ~cp_bool(row,col)
                Falpha=[Falpha, Alpha(k)==Alpha(k-1)];
            end
            k=k+1; %next alpha
        end
    end
end

% constraints for reachable set
Freach = [];
for k=1:L+1
    % x = c + G * \alpha, 
    c = reach_zonos{k}.center;
    G = reach_zonos{k}.generators;

    Freach = [Freach, X(:,k) == (c+G*Alpha(1:size(G,2))')];
end

options = Sys.solver_options;
output_controller =  {Alpha,X,Pstl};

if numel(stl_list) == 0
    Pstl = sdpvar(1,1);
end

%% Objective function, minimize robustness of stl formula
obj = sum(sum(Pstl(:,1:end))); 

milp = optimizer([Fstl, Falpha, Freach],obj,options,[], output_controller);


