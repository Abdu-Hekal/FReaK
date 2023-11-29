function Sys = setupReach(Sys)
% setupReach - Set up the evolution of the koopman linearized system 
%   using the reachable set constraints
%
% Syntax:
%    Sys = setupReach(Sys)
%
% Description:
%    This function sets up the optimization constraints for the reachable
%    set in the Koopman MILP formulation. The constraints enforce that the
%    system states follow the reachable set dynamics computed during the
%    falsification process. The function is part of the KoopMILP class and
%    is called during the setup process. Note that evolution is either described
%    directly using the setupReach function or using reachable sets (this function)
%
% Inputs:
%    Sys - KoopMILP object
%
% Outputs:
%    Sys - KoopMILP object with updated properties related to reachable
%          set optimization constraints.
%
% Example:
%    Sys = setupReach(Sys);
%
% See also: KoopMILP, koopSetupDynamics
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---


alpha=Sys.alpha; x=Sys.x; L=Sys.L;
%% Reachset constraints
Freach = [];
% Constraints for reachable set
for k=1:L+1
    % x = c + G * \alpha, 
    c = Sys.reachZonos{k}.center;
    G = Sys.reachZonos{k}.generators;
    Freach = [Freach, x(:,k) == (c+G*alpha(1:size(G,2))')];
end

%% Reachset constraints
Sys.Fdyn=Freach;

