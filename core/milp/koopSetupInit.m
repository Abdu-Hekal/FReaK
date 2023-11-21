function Sys = koopSetupInit(Sys)
% koopSetupInit - Set up the optimization constraints for the initial set
% and input constraints in the Koopman MILP formulation.
%
% Syntax:
%    Sys = koopSetupInit(Sys)
%
% Description:
%    This function sets up the optimization constraints for the initial
%    set and input constraints in the Koopman MILP formulation. It is part
%    of the KoopMILP class and is called during the setup process.
%
% Inputs:
%    Sys - KoopMILP object
%
% Outputs:
%    Sys - KoopMILP object with updated properties related to the initial
%          set and input constraints.
%
% Example:
%    Sys = koopSetupInit(Sys);
%
% See also: KoopMILP
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

