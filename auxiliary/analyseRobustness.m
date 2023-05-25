function rho = analyseRobustness(obj,x,t,offset,iterCount)
% computeRobustness - compute robustness of an STL formula on a trace
%
% Syntax:
%    rho = analyseRobustness(obj,x,t,iterCount)
%
% Inputs:
%    obj - logic formula (class stl)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%    offset- offset which is equal to robustness of trajectory on real stl
%    iterCount - trial attempt of modifying stl to reach rob==0
%
% Outputs:
%    rho - robustness of the modified formula on the trace
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

vars = getVariables(obj);
% Check input arguments
if abs(min(diff(t)) - max(diff(t))) > eps
    throw(CORAerror('CORA:notSupported',...
        'Only uniformly sampled traces are supported!'));
end

%convert stl to disjunctive normal form ??
obj = disjunctiveNormalForm(obj);
% Compute the robustness of the formula on the trace using a recursive function
global arCount
arCount=-1;
rho = recursive(obj,x,t,vars,offset,iterCount);
rho = rho(1); %robustness must return one value

end


% Recursive Function ------------------------------------------------------

function [rho] = recursive(obj,x,t,vars,offset,iterCount)
global arCount
% recursive function for computing the robustness of an STL formula on a trace
if strcmp(obj.type,'&')

    stlLhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    stlRhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    rho = min(stlLhs,stlRhs);

elseif strcmp(obj.type,'|')

    stlLhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    stlRhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    rho = max(stlLhs,stlRhs);

elseif strcmp(obj.type,'~')

    inner = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rho = -inner;

elseif strcmp(obj.type,'<')

    arCount=arCount+1;
    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rho = obj.rhs-lhs;
    if arCount == iterCount
        rho = rho - offset;
    end

elseif strcmp(obj.type,'>')

    arCount=arCount+1;
    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rho = lhs - obj.rhs;
    if arCount == iterCount
        rho = rho - offset;
    end

elseif strcmp(obj.type,'<=')

    arCount=arCount+1;
    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rho = obj.rhs-lhs - eps;
    if arCount == iterCount
        rho = rho - offset;
    end

elseif strcmp(obj.type,'>=')

    arCount=arCount+1;
    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rho = lhs - obj.rhs + eps;
    if arCount == iterCount
        rho = rho - offset;
    end

elseif strcmp(obj.type,'+')

    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    rho = lhs + rhs;

elseif strcmp(obj.type,'-')

    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    rhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    rho = lhs - rhs;

elseif strcmp(obj.type,'*')

    lhs = obj.lhs;
    rhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    rho = lhs * rhs;

elseif strcmp(obj.type,'true')

    rho = inf*ones(size(t));

elseif strcmp(obj.type,'false')

    rho = -inf*ones(size(t));

elseif strcmp(obj.type,'variable')
    var_idx = strcmp(vars,obj.var);
    rho = x(:,var_idx);

elseif strcmp(obj.type,'finally')

    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    ind = find(t >= obj.from & t <= obj.to);
    %initialize with NaN?
    rho=NaN(size(t)); cnt = 1;
    %check that time from is less than final time for traj
    if t(end) >= obj.from
        while ~isempty(ind) && ind(1) <= length(t)
            ind = ind(ind <= length(t));
            rho(cnt) = max(lhs(ind));
            cnt = cnt + 1; ind = ind + 1;
        end
    end


elseif strcmp(obj.type,'globally')
    lhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    ind = find(t >= obj.from & t <= obj.to);
    %initialize with NaN, how will this affect optim func?
    rho=NaN(size(t)); cnt = 1;
    %check that time to is less than final time for traj
    if t(end) >= obj.to
        while ~isempty(ind) && ind(end) <= length(t)
            ind = ind(ind <= length(t));
            rho(cnt) = min(lhs(ind));
            cnt = cnt + 1; ind = ind + 1;
        end
    end

elseif strcmp(obj.type,'release')

    stlLhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    stlRhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    ind = find(t >= obj.from & t <= obj.to);
    rho=NaN(size(t)); cnt = 1;

    %check that time from is less than final time for traj
    if t(end) >= obj.from
        inner = NaN(size(stlRhs));
        for k = 1:length(inner)
            inner(k) = max(stlRhs(k),max(stlLhs(1:k)));
        end
        while ~isempty(ind) && ind(1) <= length(t)
            ind = ind(ind <= length(t));
            rho(cnt) = min(inner(ind));
            cnt = cnt + 1; ind = ind + 1;
        end
    end

elseif strcmp(obj.type,'until')
    stlLhs = recursive(obj.lhs,x,t,vars,offset,iterCount);
    stlRhs = recursive(obj.rhs,x,t,vars,offset,iterCount);

    ind = find(t >= obj.from & t <= obj.to);
    rho=NaN(size(t)); cnt = 1;

    %check that time from is less than final time for traj
    if t(end) >= obj.from
        inner = NaN(size(stlRhs));
        for k = 1:length(inner)
            inner(k) = min(stlRhs(k),min(stlLhs(1:k)));
        end
        while ~isempty(ind) && ind(1) <= length(t)
            ind = ind(ind <= length(t));
            rho(cnt) = max(inner(ind));
            cnt = cnt + 1; ind = ind + 1;
        end
    end
end
%check that stl complies with times of trace
if all(isnan(rho))
    error(['Error in logic formula defined, make sure that your time' ...
        ' intervals are feasible corresponding to times of trace'])
end
end

