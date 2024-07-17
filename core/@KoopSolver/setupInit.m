function Sys = setupInit(Sys)
% setupInit - Set up the optimization constraints for the initial set
% and input constraints in the Koopman solver formulation. used only if
% direct encoding (without reachability) is used
%
% Syntax:
%    Sys = setupInit(Sys)
%
% Description:
%    This function sets up the optimization constraints for the initial
%    set and input constraints in the Koopman solver formulation. It is part
%    of the KoopSolver class and is called during the setup process.
%
% Inputs:
%    Sys - KoopSolver object
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to the initial
%          set and input constraints.
%
% Example:
%    Sys = setupInit(Sys);
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---

%setup variables 
L=Sys.L; %number of analysis points
Sys.x = sdpvar(Sys.nx+Sys.nObs, L+1); %states
Sys.u = sdpvar(Sys.nu, L+1); %inputs

U=Sys.U; u=Sys.u; 

%% Input constraints
assert(isequal(class(U),'interval'), "U must be an interval if not using reachability (in the version)")
F= repmat(U.inf,1,L+1)<=u<=repmat(U.sup,1,L+1);

Sys.Finit=F;

