function Sys = addPredConstr(Sys,predTimeConstrs,nPreds,hardcoded)
% addPredConstr - add constraints on predicates of an STL formula at
% desired times.
%
% Syntax:
%    Sys = addPredConstr(Sys, hardcoded)
%
% Description:
%    This function adds (soft) constraints for the predicates in Signal
%    Temporal Logic (STL) formula in the Koopman Solver formulation. The
%    constraints ensure that the predicate is satisfied at the desired time
%    point, the fcn also adds a variable Ostl which is set as objective fcn
%    to maximize satisfaction of predicates (or minimize dissatisfaction),
%    this variables is why we define the constraints as soft constraint.
%    e.g. for a predicate x>5, we constrain x>5+b, and set the optimization
%    problem to maximize b.
%
% Inputs:
%    Sys - KoopSolver object
%    predTimeConstrs - a dictionary where keys are indices of predicate to
%                      add constraint to and values are time points where 
%                      predicates should be constrained.
%    nPreds - number of predicates in the stl formula
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
%    Sys = addPredConstr(Sys, true);
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     10-March-2024
% Last update: ---
% Last revision: ---

%if first time adding constraints, setup robustness param and offset params
if isempty(Sys.Ostl)
    Sys.Pstl=sdpvar(1,1);
    if ~hardcoded
        Sys.Ostl = sdpvar(1,nPreds);
    end
end

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; 

end


