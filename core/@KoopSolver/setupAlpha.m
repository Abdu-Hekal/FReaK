function Sys = setupAlpha(Sys)
% setupAlpha - Set up the optimization variables and constraints for the
% alpha parameters in the Koopman solver formulation.
%
% Syntax:
%    Sys = setupAlpha(Sys)
%
% Description:
%    This function sets up the optimization variables and constraints for
%    the alpha parameters in the Koopman solver formulation. The function is
%    part of the KoopSolver class and is called during the setup process.
%
% Inputs:
%    Sys - KoopSolver object
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to alpha
%          optimization variables and constraints.
%
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---


%% setup
L=Sys.L; % horizon (# of steps)
reachZonos=Sys.reachZonos;

%% System dimensions and variables
nx=Sys.nx; %number of states
% variables
Sys.x = sdpvar(nx, L+1); %states
alpha = sdpvar(1, size(reachZonos{end}.generators,2));
if ~isempty(Sys.U)
    alphaU = alpha(size(reachZonos{1}.generators,2)+1:end);
    U = zonotope(Sys.U); c_u = center(U); G_u = generators(U);
    alphaU = reshape(alphaU,[size(G_u,2),length(alphaU)/size(G_u,2)]);
    c_u_ = repmat(c_u,1,size(alphaU,2));

    %append empty sdpvar for consistent length with states X
    Sys.u = [c_u_ + G_u*alphaU, sdpvar(size(c_u,1),size(c_u,2))];
end

%constraints for alpha
Falpha= -1<=alpha<=1;
%assign optim variables and outputs to system
Sys.Finit=Falpha; Sys.alpha=alpha;




