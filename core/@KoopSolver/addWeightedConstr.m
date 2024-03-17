function Sys = addWeightedConstr(Sys,predTimeConstrs,preds,hardcoded,offsetStrat)
% addWeightedConstr - add constraints on predicates of an STL formula at
% desired times.
%
% Syntax:
%    Sys = addWeightedConstr(Sys, hardcoded)
%
% Description:
%    This function adds (soft) weighted constraints for the predicates in Signal
%    Temporal Logic (STL) formula in the Koopman Solver formulation. The
%    constraints first aim to maximize satisfaction of all predicates, with
%    different weights depending on the logic operator. The weights can be
%    manipulated based on information from predTimeConstrs that records
%    critical times and corresponding predicates, where formula is least
%    (dis)satisfied
%
% Inputs:
%    Sys - KoopSolver object
%    predTimeConstrs - a dictionary where keys are indices of predicate to
%                      add constraint to and values are time points where 
%                      predicates should be constrained.
%    preds - list of all predicates
%    hardcoded - Boolean flag indicating whether the STL formula is
%                hardcoded (true) or not (false). if false, an optimizer
%                object is used.
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to the STL
%          formula optimization constraints, Namely:
%     Fstl: Yalmip constraints on predicates
%     Pstl:  a struct containing YALMIP decision variables representing
%           the quantitative satisfaction of the predicate(s)
%     Ostl:  a struct containing YALMIP parameter variables representing the
%       offset for each inequality.
%
% Example:
%    Sys = addWeightedConstr(Sys, true);
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     10-March-2024
% Last update: ---
% Last revision: ---

%if first time adding constraints, setup robustness param and offset params
if isempty(Sys.Pstl)
    Sys.Pstl=sdpvar(1,1);
end
if isempty(Sys.Ostl) && ~hardcoded
    Sys.Ostl = sdpvar(1,numel(preds));
end

%setup variables
x=Sys.x; u=Sys.u;
pstl=Sys.Pstl; ostl=Sys.Ostl;

Sys.solverTimePoints

error('stop')

end


