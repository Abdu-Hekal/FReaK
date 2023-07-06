function [idx,sign,mus]=bReachCulprit(Bdata,phi,rob)
%function that gets idx of predicate (subformula) that is responsible for robustness value

mus = STL_ExtractPredicates(phi);
% Obtain unique values and their counts
uniqueMus={};
counts = {};
idx=-1;
sign=0;
signs=[1,-1];
% no point of offset if only one predicate
if numel(mus) <= 1
    return
end
for ii=1:numel(mus)
    for jj=1:numel(signs)
        pred = mus(ii);
        predStl=regexprep(evalc('display(pred)'), '\n', '');
        predStl = replace(predStl,'(','');
        predStl = replace(predStl,')','');
        modPredStl = strcat(predStl,'+',char(vpa(signs(jj)*rob)));
        stl=regexprep(evalc('display(phi)'), '\n', '');
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
        Rphi = BreachRequirement(modPhi);
        newRob=Rphi.Eval(Bdata);
        if newRob<rob
            idx = ii;
            sign=signs(jj);
            return
        end
    end
end
end
