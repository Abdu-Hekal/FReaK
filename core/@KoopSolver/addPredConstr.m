function Sys = addPredConstr(Sys,predTimeConstrs,preds,hardcoded,offsetStrat)
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
%    Sys = addPredConstr(Sys, true);
%
% See also: KoopSolver
%
% Author:      Abdelrahman Hekal
% Written:     10-March-2024
% Last update: ---
% Last revision: ---

%The following 3 lines ensure only unique constraints exist,
% i.e. only new predicate constraints are added
% Use cellfun with anonymous functions to convert structs to strings
strCell = cellfun(@(s) jsonencode(s), predTimeConstrs, 'UniformOutput', false);
uniqueStrCell = unique(strCell, 'stable');
% Convert the unique cell array of strings back to a cell array of structs
predTimeConstrs = cellfun(@jsondecode, uniqueStrCell, 'UniformOutput', false);

%if first time adding constraints, setup offset params
if isempty(Sys.Ostl) && ~hardcoded
    Sys.Ostl = sdpvar(1,numel(preds));
end

%setup variables
x=Sys.x; u=Sys.u;
pstl=Sys.Pstl; ostl=Sys.Ostl;

%if hardcoded (no optimizer object) and we use offset, then we need to encode all constrs from scratch
if hardcoded && abs(offsetStrat)
    constrs=predTimeConstrs;
    Sys.Fstl=[]; %empty current constraints
else %only add new constrs
    constrs=predTimeConstrs(length(Sys.Fstl)+1:end);
end

for p=1:numel(preds)
    pred=preds{p};
    %for each predicate, find all times where constraint must hold
    matchingIndices = cellfun(@(x) strcmp(x.pred,pred), constrs);
    %we add 1, as time indices start from 1 instead of 0
    allTimes = cellfun(@(x) x.time+1, constrs(matchingIndices));
    if ~isempty(allTimes)
        %remove new line at end of pred if any
        pred=regexprep(pred, '\n', '');
        %remove enclosing brackets of pred
        pred = regexprep(pred, '^(\()|\)$', '');

        pred = regexprep(pred, '\[t\]', ''); %remove [t]
        pred = regexprep(pred,'x(\w*)','x($1,t)'); %add time to state variables
        pred = regexprep(pred,'u(\w*)','u($1,t)'); %add time to inputs
        pred = regexprep(pred, '\<t\>', sprintf('[%s]', strjoin(arrayfun(@num2str, unique(allTimes), 'UniformOutput', false), ', ')));

        %replace all with strict inequality, hack to later remove all
        pred = replace(pred,'<=','<');
        pred = replace(pred,'>=','>');

        %flip signs, as we constrain satisfaction not falsification
        pred = regexprep(pred, '>', 'TEMP_STRING');
        pred = regexprep(pred, '<', '>');
        pred = regexprep(pred, 'TEMP_STRING', '<');

        if hardcoded
            if numEntries(Sys.offsetMap)>0 && isKey(Sys.offsetMap,p)
                offset=Sys.offsetMap(p);
            else
                offset=0;
            end
            if contains(pred, '<')
                pred=strcat(pred,' + pstl - offset');
            elseif contains(pred, '>')
                pred=strcat(pred,' - pstl - offset');                          
            end
        else
            if contains(pred, '<')
                pred=strcat(pred,' + pstl - ostl(p)');
            elseif contains(pred, '>')
                pred=strcat(pred,' - pstl - ostl(p)');
            end
        end

        %remove any strict inequalities
        pred = replace(pred,'<','<=');
        pred = replace(pred,'>','>=');
        
        z_eval = eval(pred);
        Sys.Fstl=[Sys.Fstl,z_eval];
    end
end
end


