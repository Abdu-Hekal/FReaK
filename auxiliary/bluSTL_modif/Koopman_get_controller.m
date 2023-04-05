function controller = Koopman_get_controller(Sys,~)
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
ts=Sys.ts; % sampling time
L=Sys.L;   % horizon (# of steps)

%% System dimensions and variables
nu=Sys.nu;
nx=Sys.nx;

% variables
X = sdpvar(nx, L+1);
U = sdpvar(nu, L);

%% STL formula
Fstl = [];
Pstl = [];
varStd = struct('X',X,'U',U);

if isstruct(Sys.var)
    %remove overlapping fields from std
    var = rmfield(varStd, intersect(fieldnames(Sys.var), fieldnames(varStd)));
    keys = [fieldnames(var); fieldnames(Sys.var)];
    var = cell2struct([struct2cell(varStd); struct2cell(Sys.var)], keys, 1);
else
    var = varStd;
end

stl_list= STLC_parse_stl_labels(Sys);
M = Sys.bigM;


for i = 1:numel(stl_list)
    phi = STLformula('phi', stl_list{i});

    %what time steps to compute robustness for? does this work? why?
    [Fphi, Pphi] = Koopman_MILP_robust(phi, 1, L+1, ts, var,M);
    
    Pstl = [Pstl; Pphi];
    Fstl = [Fstl Fphi];

end


%% Input constraints
Fu = [];

% Bounds
for iu = 1:nu
    Fu = [ Fu, Sys.u_lb(iu) <= U(iu,:) <= Sys.u_ub(iu)] ;  % bounds constraints on u
end


%% Dynamics constraints
Fdyn = [];

[Ad,Bd]=ssdata(Sys.sysd);

Bdu=Bd(:,1:nu);

% Constraints for states (if any)
for k=1:L+1
    if k==1
        if size(Sys.x0,2) == 1
            Fdyn = [Fdyn, Sys.x0(:,1) <= X(:,k) <= Sys.x0(:,1)];  % bounds constraints on initial state
        else
            Fdyn = [Fdyn, Sys.x0(:,1) <= X(:,k) <= Sys.x0(:,2)];  % bounds constraints on initial state
        end

    else
        %dynamics
        Fdyn = [Fdyn, X(:,k) == (Ad*X(:,k-1) + Bdu*U(:,k-1))];
    end
end

options = Sys.solver_options;
output_controller =  {U,X,Pstl};


if numel(stl_list) == 0
    Pstl = sdpvar(1,1);
end

%% Objective function
obj = sum(sum(Pstl(:,1:end))); 

controller = optimizer([Fdyn, Fstl, Fu],obj,options,[], output_controller);


