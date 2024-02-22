function Sys = setupStl(Sys,hardcoded)
% setupStl - Set up the optimization constraints for the robustness of
%   the STL formula in the Koopman MILP formulation.
%
% Syntax:
%    Sys = setupStl(Sys, hardcoded)
%
% Description:
%    This function sets up the optimization constraints for the Signal
%    Temporal Logic (STL) formula in the Koopman MILP formulation. The
%    constraints ensure that the robustness of the STL formula is
%    minimized during the falsification process. The function is part of
%    the KoopMILP class and is called during the setup process.
%
% Inputs:
%    Sys - KoopMILP object
%    hardcoded - Boolean flag indicating whether the STL formula is
%                hardcoded (true) or not (false). if false, an optimizer
%                object is used.
%
% Outputs:
%    Sys - KoopMILP object with updated properties related to the STL
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
% See also: KoopMILP
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
%evaluate stl formula at specified time indices, note +1 is used to start
%from 1 instead of 0
timeIdxs = floor(Sys.solverTimePoints/Sys.koopdt)+1;
x = Sys.x(:,timeIdxs);
u = Sys.u(:,timeIdxs);
var = struct('x',x,'u',u);
L=size(x,2);

phi= Sys.stl;
M = Sys.bigM;

if Sys.normalize
    normz = Sys.normz;
else
    normz = [];
end

if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=0;
end

[Fstl, Pstl, Ostl] = koopStl(phi,1,L,Sys.solverTimePoints,var,M,normz,hardcoded,Sys.offsetMap);

%assign stl optim variables and constraints
Sys.Fstl=Fstl; Sys.Pstl=Pstl; Sys.Ostl=Ostl;

end


