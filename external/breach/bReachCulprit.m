function [critPreds, critTimes, preds]=bReachCulprit(Bdata,set,sepCaluses)
% bReachCulprit - Finds indices of critical predicates and critical time
% points
%
% Note that unlike typical critical pred computation where the goal is to
% find the predicate responsible for the current robustness value, the
% functions finds all predicates that are responsible for a +ve value for
% this stl and trajectory, (i.e. not falsified)
%
% Syntax:
%   [critPreds, critTimes, preds] = bReachCulprit(Bdata, set)
%
% Description:
%   This function iterates through the clauses of the provided STL formula set
%   and identifies predicates responsible for positive robustness values.
%   It also can return corresponding critical time points, i.e. where the
%   trajectory violates the STL the most
%
% Inputs:
%   Bdata - Data for the Breach system trajectory.
%   set - spec set containing STL formula.
%   sepClauses - boolean to select whether to seperate clauses or not.
%
% Outputs:
%   critPreds - dictionary of indices of predicates responsible for positive
%               robustness values. Key: predicate index, Value: offset.
%   critTimes - cell array of structs which store critical predicates and their
%               corresponding critical times
%   preds -     list of predicates as strings
%
% See also: getClauses, coraBreachConvert, recursiveOffset
%
% Author: Abdelrahman Hekal
% Written: November 2023
% Last update: 13-March-2024
% Last revision: ---

critPreds = dictionary();
critTimes = {};
preds = {};

idx=0;
if nargin<=2 || sepCaluses
    clauses = getClauses(set);
else
    clauses={set};
end
for ij=1:numel(clauses)
    clause = clauses{ij};
    stl=coraBreachConvert(clause);
    phi = STL_Formula('phi',stl);

    [critPreds, critTimes] = recursiveOffset(critPreds,critTimes,idx,phi,Bdata);
    mus = STL_ExtractPredicates(phi);
    for ix=1:numel(mus)
        idx = idx+1; %increase idx, i.e. one more predicate
        % add predicates to map of predicates
        preds{idx} = get_st(mus(ix));
    end
end
% no point of offset if only one clause with one predicate, return empty
% critPreds
if numel(clauses) <= 1
    if numel(mus) <= 1
        critPreds = dictionary();
        return
    end
end
end

function [critPreds, critTimes] = recursiveOffset(critPreds,critTimes,idx,phi,Bdata)
%compute current robustness and extract all predicates
traj=Bdata.P.traj{1};
P=Sselect(Bdata.P,1);
val = STL_Eval(Bdata.Sys, phi, P, traj, traj.time);
rob=val(1);

mus = STL_ExtractPredicates(phi);
% Obtain unique values and their counts
uniqueMus={};
counts = {};
for ii=1:numel(mus)
    if numEntries(critPreds)>0 && isKey(critPreds, idx+ii)
        continue; % Skip if pred is already offset
    end
    pred = mus(ii);
    %remove new line at end of pred if any
    predStl=regexprep(get_st(pred), '\n', '');
    %remove enclosing brackets of pred
    predStl = regexprep(predStl, '^(\()|\)$', '');
    %if in cnf, no need to iterate over signs, sign is +ve is pred>c 
    % or -ve if pred<c, if there exist ~, might need
    %to iterate oversigns or handle negation.
    if contains(predStl, '>')
        sign=1;
    elseif contains(predStl, '<')
        sign=-1;
    end
    %modify predicate with current value of rob
    modPredStl = strcat(predStl,'+',char(vpa(sign*rob)));
    stl=regexprep(disp(phi), '\n', '');
    pattern = strcat('(.*?)',regexptranslate('escape', predStl),'(.*?)');
    %only count once for each pred
    index = find(strcmp(predStl, uniqueMus));
    if isempty(index)
        uniqueMus{end+1} = predStl;
        counts{end+1} = 1;
        count=1;
    else
        counts{index} = counts{index}+1;
        count=counts{index};
    end
    modStl = regexprep(stl,pattern,strcat('$1',modPredStl,'$2'),count);
    modPhi = STL_Formula('phi',modStl);

    val = STL_Eval(Bdata.Sys, modPhi, P, traj, traj.time);
    newRob=val(1);

    if newRob<rob
        %get critical time point by computing robustness of predicate
        %at each time point. This also checks that the pred found is in
        %fact the critical pred.
        val = STL_Eval(Bdata.Sys, pred, P, traj, traj.time);
        timeIdx=find(abs(val-rob)<1e-13,1); %val==rob with tolerance of 1e-13
        if ~isempty(timeIdx)
            critPreds(idx+ii) = sign*rob;
            %store crit time struct
            critTime=struct;
            critTime.time=traj.time(timeIdx);
            critTime.pred=get_st(pred);
            critTime.weight=rob;
            critTimes{end+1}=critTime;
            if newRob > 1e-13 %not yet offset all responsible predicates
                [critPreds,critTimes] = recursiveOffset(critPreds,critTimes,idx,modPhi,Bdata);
            else %found all responsible predicates
                return
            end
        end
    end
end
end

