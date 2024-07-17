function Sys = setupCP(Sys)
% setupCP - Set up constraints on control points for the Koopman solver (KoopSolver) object.
%
% Syntax:
%    Sys = setupCP(Sys)
%
% Description:
%    This function sets up constraints on control points for the Koopman
%    solver (KoopSolver) object. The constraints enforce that the input is
%    activated only at the control points specified by the boolean array
%    Sys.cpBool. The remaining input points (if exist) are set to be equal 
%    to previous control point (piecewise-constant-interpolation). If
%    Sys.cpBool is empty or all true, no additional constraints are added.
%
% Inputs:
%    Sys - Koopman solver (KoopSolver) object
%
% Outputs:
%    Sys - Updated Koopman solver (KoopSolver) object with control point constraints set up
%
% Example:
%    Sys = setupCP(Sys);
%
% See also: optimize
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
%------------- BEGIN CODE --------------

% TODO: constrain to same value as prev only makes sense
% for piecewise-constant inputs. can we consider other interpolation
% methods?
%% constraints on input based on boolean of control points
cpBool=Sys.cpBool; %boolean of control points
if ~isempty(cpBool) && ~all(cpBool,'all') %piecewise constant signal and not pulse. if all cpbool is ones, then no constraints needed.
    F_cp = Sys.u(:,1:end-1) .* ~cpBool' == [zeros(size(Sys.u,1),1),Sys.u(:,1:end-2)] .* ~cpBool';
    Sys.Finit=[Sys.Finit, F_cp];
end
end