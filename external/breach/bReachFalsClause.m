function [culpritSet,nClauses]=bReachFalsClause(Bdata,set)
%function that get subclause with least robustness value from an stl in
%conjunctive normal form

clauses = getClauses(set);
nClauses = numel(clauses);
rob=inf;
%only one clause
if numel(clauses) <=1
    culpritSet=set;
    return;
end
for ij=1:numel(clauses)
    clause = clauses{ij};
    stl=replace(coraBlustlConvert(clause),"(t)","[t]");
    phi = STL_Formula('phi',stl);

    %compute current robustness
    Rphi = BreachRequirement(phi);
    rob_=Rphi.Eval(Bdata);

    if rob_<rob
        culpritSet = clause;
        rob=rob_;
    end
end
end
