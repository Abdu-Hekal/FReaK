function [critPreds, critTimes]=bReachCulprit(Bdata,set)
% bReachCulprit - Finds indices of critical predicates and critical time
% points
%
% Note that unlike typical critical pred computation where the goal is to
% find the predicate responsible for the current robustness value, the
% functions finds all predicates that are responsible for a +ve value for
% this stl and trajectory, (i.e. not falsified)
%
% Syntax:
%   [critPreds] = bReachCulprit(Bdata, set)
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
%
% Outputs:
%   critPreds - Map of indices of predicates responsible for positive
%               robustness values. Key: predicate index, Value: offset.
%
% Example:
%   [critPreds] = bReachCulprit(Bdata, set);
%
% See also: getClauses, coraBreachConvert, recursiveOffset
%
% Author: Abdelrahman Hekal
% Written: [Date]
% Last update: [Date]
% Last revision: ---

critPreds = containers.Map('KeyType', 'double', 'ValueType', 'double');
critTimes = containers.Map('KeyType', 'double', 'ValueType', 'double');
idx=0;
clauses = getClauses(set);
for ij=1:numel(clauses)
    clause = clauses{ij};
    stl=coraBreachConvert(clause);
    phi = STL_Formula('phi',stl);

    [critPreds, critTimes] = recursiveOffset(critPreds,critTimes,idx,phi,Bdata);
    mus = STL_ExtractPredicates(phi);
    idx = idx+numel(mus); %increase idx to search for next mus
end
% no point of offset if only one clause with one predicate, return empty
% critPreds
if numel(clauses) <= 1
    if numel(mus) <= 1
        critPreds = containers.Map('KeyType', 'double', 'ValueType', 'double');
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
signs=[1,-1];
for ii=1:numel(mus)
    if isKey(critPreds, idx+ii)
        continue; % Skip if pred is already offset
    end
    pred = mus(ii);
    %remove new line at end of pred if any
    predStl=regexprep(disp(pred), '\n', '');
    %remove enclosing brackets of pred
    predStl = regexprep(predStl, '^(\()|\)$', '');
    %iterate over signs. TODO: if in cnf, no need to iterate over signs,
    %sign is +ve is pred>c or -ve if pred<c
    for jj=1:numel(signs)
        %modify predicate with current value of rob
        modPredStl = strcat(predStl,'+',char(vpa(signs(jj)*rob)));
        stl=regexprep(disp(phi), '\n', '');
        pattern = strcat('(.*?)',regexptranslate('escape', predStl),'(.*?)');
        %only count once for each pred
        if jj==1
            index = find(strcmp(predStl, uniqueMus));
            if isempty(index)
                uniqueMus{end+1} = predStl;
                counts{end+1} = 1;
                count=1;
            else
                counts{index} = counts{index}+1;
                count=counts{index};
            end
        end
        modStl = regexprep(stl,pattern,strcat('$1',modPredStl,'$2'),count);
        modPhi = STL_Formula('phi',modStl);

        val = STL_Eval(Bdata.Sys, modPhi, P, traj, traj.time);
        newRob=val(1);

        if newRob<rob
            critPreds(idx+ii) = signs(jj)*rob;
            %get critical time point by computing robustness of predicate
            %at each time point
            val = STL_Eval(Bdata.Sys, pred, P, traj, traj.time);
            timeIdx=find(abs(val-rob)<1e-13,1); %val==rob with tolerance of 1e-13
            if ~isempty(timeIdx)
                critTimes(idx+ii) = traj.time(timeIdx);
            end
            if newRob > 1e-13 %not yet offset all responsible predicates
                [critPreds] = recursiveOffset(critPreds,critTimes,idx,modPhi,Bdata);
            else %found all responsible predicates
                return
            end
        end
    end
end
end

