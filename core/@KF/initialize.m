function [obj,trainset,soln,specSolns] = initialize(obj)
% INITIALIZE Initialize the Koopman Falsification (KF) object and check for required parameters.
%
% Syntax:
%    [obj, trainset] = initialize(obj)
%
% Description:
%    This function initializes the Koopman Falsification (KF) object by checking
%    and setting the required parameters. It also ensures that autokoopman is
%    installed and imported in the Python environment. The function returns the
%    initialized object and an empty training set structure.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%    obj       - Initialized Koopman Falsification (KF) object
%    trainset  - Empty training set structure
%    soln      - struct to store solution
%    specSolns - dictionary to store previous solution for each spec
%
% Example:
%    [obj, trainset] = initialize(obj);
%
% See also: falsify
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
% -------------------------- Auxiliary Functions --------------------------
% SETCPBOOL Set control point boolean array for the Koopman Falsification (KF) object.
%
% Syntax:
%    obj = setCpBool(obj)
%
% Description:
%    This function sets the control point boolean array for the Koopman Falsification (KF) object.
%    The boolean array represents the control points in the Koopman model.
%
% Inputs:
%    obj - Koopman Falsification (KF) object
%
% Outputs:
%    obj - Updated Koopman Falsification (KF) object with the control point boolean array
%
% Example:
%    obj = setCpBool(obj);
%
% See also: initialize
%
% Author:      Abdelrahman Hekal
% Written:     19-November-2023
% Last update: ---
% Last revision: ---
% -----------------------------------------------------------------------------
% ------------- BEGIN CODE --------------

%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

assert(isa(obj.model, 'string') | isa(obj.model,"char")| isa(obj.model,'function_handle'), 'obj.model must be a (1)string: name of simulink obj or a (2)function handle')
assert(isa(obj.R0, 'interval'), 'Initial set (obj.R0) must be defined as an CORA interval')
assert(~isempty(obj.T) & isnumeric(obj.T), 'Time horizon (obj.T) must be defined as a numeric')
assert(~isempty(obj.dt) & isnumeric(obj.dt), 'Time step (obj.dt) must be defined as a numeric')
assert(isa(obj.spec, 'specification'), 'Falsifying spec (obj.spec) must be defined as a CORA specification')
all_steps = obj.T/obj.dt;
assert(floor(all_steps)==all_steps,'Time step (dt) must be a factor of Time horizon (T)')

%set autokoopman timestep if it is not set, else check it is compliant.
if ~isfield(obj.ak,'dt')
    obj.ak.dt=obj.dt;
else
    allAbstrSteps = obj.T/obj.ak.dt;
    assert(floor(allAbstrSteps)==allAbstrSteps,'Time step of koopman (ak.dt) must be a factor of Time horizon (T)')
end

%set solver timestep if it is not set, else check it is compliant.
if ~isfield(obj.solver,'dt')
    obj.solver.dt=obj.ak.dt;
else
    abstr = obj.solver.dt/obj.ak.dt; %define abstraction ratio
    assert(floor(abstr)==abstr,'Time step of solver (solver.dt) must be a multiple of koopman time step (ak.dt)')
    allAbstrSteps = obj.T/obj.solver.dt;
    assert(floor(allAbstrSteps)==allAbstrSteps,'Time step of solver (solver.dt) must be a factor of Time horizon (T)')
end

% ensure that autokoopman rank is an integer
obj.ak.rank=int64(obj.ak.rank);

% clear yalmip
yalmip('clear')

if ~isempty(obj.U) %check if obj has inputs
    assert(isa(obj.U, 'interval'), 'Input (obj.U) must be defined as an CORA interval')
    %if no control points defined, set as control point at every step ak.dt
    if isempty(obj.cp)
        obj.cp=obj.T/obj.ak.dt*ones(1,length(obj.U));
    end

    assert(length(obj.U)==length(obj.cp),'Number of control points (obj.cp) must be equal to number of inputs (obj.U)')
end

%set cpBool
obj=setCpBool(obj);

%reset struct to store soln
soln=struct;
soln.falsified=false;
soln.koopTime=0; soln.reachTime=0;
soln.optimTime=0; soln.simTime=0;
soln.sims=0;
soln.best.rob=inf;
%reset dict to store prev soln for each spec
specSolns = dictionary(obj.spec,struct);
%empty struct to store training data
trainset.X = {}; trainset.XU={}; trainset.t = {};

end

% -------------------------- Auxiliary Functions --------------------------

function obj=setCpBool(obj)
all_steps = obj.T/obj.ak.dt;
obj.cpBool = zeros(all_steps,length(obj.U));
for k=1:length(obj.cp)
    if ~isempty(obj.U)
        assert(obj.cp(k)>0, 'your model has inputs defined, number of control points must be greater than zero')
    end
    if obj.cp(k) == 1
       assert(obj.cp(k)>0, 'if number of control points is 1, pconst interpolation must be used')
    end
    step = (obj.T/obj.ak.dt)/obj.cp(k);
    step = max(1,step); %if step<1, then we have more control points than ak steps, set control points equal to number of steps
    assert(floor(step)==step,'number of control points (cp) must be a factor of T/ak.dt')
    obj.cpBool(1:step:end,k) = 1;
    if ~isempty(obj.cpBool) && ~all(obj.cpBool(:,k)) %if cpbool is not just ones, then interpolation scheme must be pconst
        assert(strcmp(obj.inputInterpolation,'previous'),'Currently only an input interpolation of "previous" is supported for a number of control points less than T/ak.dt')
    end
end
end