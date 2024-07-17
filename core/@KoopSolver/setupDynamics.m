function Sys = setupDynamics(Sys)
% setupDynamics - Set up the optimization constraints for the system
% dynamics in the Koopman solver formulation.
%
% Syntax:
%    Sys = setupDynamics(Sys)
%
% Description:
%    This function sets up the optimization constraints for the system
%    dynamics in the Koopman solver formulation. The dynamics constraints
%    enforce the linear evolution of the system states over the prediction
%    horizon. The function is part of the KoopSolver class and is called
%    during the setup process. Note that evolution is either described
%    directly using this function or using reachable sets (see koopSetupReach)
%
% Inputs:
%    Sys - KoopSolver object
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to dynamics
%          optimization constraints.
%
% Example:
%    Sys = setupDynamics(Sys);
%
% See also: KoopSolver, koopSetupReach
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---


L=Sys.L;

A=Sys.A; B=Sys.B;
x=Sys.x; u=Sys.u; 


%% Dynamics constraints
assert(all(Sys.X0.inf==Sys.X0.sup),'Initial region "R0" must be a point if reachability is not used')
F = x(:,1) == Sys.g(center(Sys.X0));
for k=2:L+1
    % x = Ax + Bu 
    F = [F, x(:,k) == A*x(:,k-1) + B*u(:,k-1)];
end

%% Dynamics constraints
Sys.Fdyn=F;

