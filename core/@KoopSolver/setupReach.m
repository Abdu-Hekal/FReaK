function Sys = setupReach(Sys)
% setupReach - Set up the evolution of the koopman linearized system 
%   using the reachable set constraints
%
% Syntax:
%    Sys = setupReach(Sys)
%
% Description:
%    This function sets up the optimization constraints for the reachable
%    set in the Koopman solver formulation. The constraints enforce that the
%    system states follow the reachable set dynamics computed during the
%    falsification process. The function is part of the KoopSolver class and
%    is called during the setup process. Note that evolution is either described
%    directly using the setupReach function or using reachable sets (this function)
%mi
% Inputs:
%    Sys - KoopSolver object
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to reachable
%          set optimization constraints.
%
% Example:
%    Sys = setupReach(Sys);
%
% See also: KoopSolver, koopSetupDynamics
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---


alpha=Sys.alpha; x=Sys.x;
%% Reachset constraints
Freach = [];
%set constraint on x every n koopman time steps. Note, no need for other
%constraints on x, because stl is evaluated on x every n steps, subject to
%solver dt
timeIdxs = floor(Sys.solverTimePoints/Sys.koopdt)+1; 
% Constraints for reachable set
for k=timeIdxs
    % x = c + G * \alpha, 
    c = Sys.reachZonos{k}.center;
    G = Sys.reachZonos{k}.generators;
    Freach = [Freach, x(:,k) == (c+G*alpha(1:size(G,2))')];
end

%% Reachset constraints
Sys.Fdyn=Freach;

