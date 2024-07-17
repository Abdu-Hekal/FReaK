function Sys = setupStl(Sys,hardcoded,weighted)
% setupStl - Set up the optimization constraints for the robustness of
%   the STL formula in the Koopman Solver formulation.
%
% Syntax:
%    Sys = setupStl(Sys, hardcoded)
%
% Description:
%    This function sets up the optimization constraints for the Signal
%    Temporal Logic (STL) formula in the Koopman Solver formulation. The
%    constraints are either milp constraints that minimize robustness of 
%    stl formula or weighted predicate constraints (LP), where the weights
%    are iteratively modified in the falsificaiton framework to favour the
%    satisfaction of critical preds at crit times moments. If predTimeConstrs
%    and preds are passed, then the latter approach is used, else milp
%    approach is used. This function is part of the KoopSolver class and is 
%    called during the setup process.
%
% Inputs:
%    Sys - KoopSolver object
%    hardcoded - Boolean flag indicating whether the STL formula is
%                hardcoded (true) or not (false). if false, an optimizer
%                object is used.
%    predTimeConstrs (optional) - a dictionary where keys are indices of predicate to
%                      add constraint to and values are time points where 
%                      predicates should be constrained.
%    preds (optional) - list of all predicates
%
% Outputs:
%    Sys - KoopSolver object with updated properties related to the STL
%          formula optimization constraints, Namely:
%     Fstl: Yalmip constraints to compute recusrive robustness of stl
%     Pstl:  a struct containing YALMIP decision variables representing
%           the quantitative satisfaction of phi over each time step in
%           kList
%     Ostl:  a struct containing YALMIP parameter variables representing the
%       offset for each inequality.
%
% Example:
%    Sys = setupStl(Sys, true);
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%if solverTimePoints is empty, then no objective fcn of stl robustness, and prob is solver
%only input/dynamics constraints
if isempty(Sys.solverTimePoints)
    Sys.Pstl=sdpvar(1,1); %dummy empty sdpvar for returned soln
    return
end


%% STL formula
phi= Sys.stl;
M = Sys.bigM;
%evaluate stl formula at specified time indices, note +1 is used to start
%from 1 instead of 0
timeIdxs = floor(Sys.solverTimePoints/Sys.koopdt)+1;
x = Sys.x(:,timeIdxs);
%check if there exists input
if ~isempty(Sys.u)
    u = Sys.u(:,timeIdxs);
else
    u=[];
end
var = struct('x',x,'u',u);
L=size(x,2);

if Sys.normalize
    normz = Sys.normz;
else
    normz = [];
end

if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
end

if nargin<=2 || ~weighted
    [Fstl, Pstl, Ostl] = koopMilpStl(phi,1,L,Sys.solverTimePoints,var,M,normz,hardcoded,Sys.offsetMap);
else
    %select weights at specified solver time points 
    weights=Sys.weights(:,timeIdxs);
    %make sure weights aren't too big
    weights=min(weights,Sys.maxWeight);
    [Fstl, Pstl, Ostl, Wstl] = koopWeightStl(phi,1,L,Sys.solverTimePoints,var,normz,hardcoded,Sys.offsetMap,weights);
    Sys.Wstl=Wstl; %set weights paramaters for optimizer object
end

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;

end


