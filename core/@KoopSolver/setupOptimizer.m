function Sys = setupOptimizer(Sys,options)
% setupOptimizer - Set up the optimizer for the Koopman solver (KoopSolver) object.
%
% Syntax:
%    Sys = setupOptimizer(Sys, options)
%
% Description:
%    This function sets up the optimizer for the Koopman solver (KoopSolver) object.
%    It defines the optimization problem, including constraints, objective,
%    and output variables. The optimizer object is stored in the Sys object
%    for later use in the optimization process.
%
% Inputs:
%    Sys - Koopman solver (KoopSolver) object
%    options - Options for the optimization process
%
% Outputs:
%    Sys - Updated Koopman solver (KoopSolver) object with the optimizer set up
%
% Example:
%    Sys = setupOptimizer(Sys, options);
%
% See also: optimize
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

constraints=[Sys.Finit, Sys.Fstl, Sys.Fdyn];
objective = Sys.Pstl; %objective is to minimize robustness of stl formula (falsification)
if ~isempty(Sys.reachZonos)
    if ~isempty(Sys.u)
        output = {Sys.x,Sys.Pstl,Sys.alpha,Sys.u};
    else
        output = {Sys.x,Sys.Pstl,Sys.alpha};
    end
else
    output = {Sys.x,Sys.Pstl,Sys.u};
end
if isempty(Sys.Wstl)
    params=Sys.Ostl;
else
   params={Sys.Ostl,Sys.Wstl};
end
% setup optimizer
Sys.optimizer = optimizer(constraints,objective,options,params,output);
end