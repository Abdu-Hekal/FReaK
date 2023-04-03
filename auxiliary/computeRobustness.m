function rho = computeRobustness(obj,x,t)
% computeRobustness - compute robustness of an STL formula on a trace
%
% Syntax:
%    rho = computeRobustness(obj,x,t)
%
% Inputs:
%    obj - logic formula (class stl)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%
% Outputs:
%    rho - robustness of the formula on the trace
%
% Example:
%    x = stl('x',2)
%    eq = globally(x(2) < -0.5,interval(0,1));
%
%    phi = -pi/2:0.01:0;
%    x = [cos(phi'),sin(phi')];
%    t = linspace(0,1,length(phi))';
%
%    rho = computeRobustness(eq,x,t)
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

% Compute the robustness of the formula on the trace using a recursive function
rho = recursive(obj,x,t,vars);
rho = rho(1); %robustness must return one value

end


% Recursive Function ------------------------------------------------------

function rho = recursive(obj,x,t,vars)
% recursive function for computing the robustness of an STL formula on a trace

if strcmp(obj.type,'&')

    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    rho = min(lhs,rhs);

elseif strcmp(obj.type,'|')

    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    rho = max(lhs,rhs);

elseif strcmp(obj.type,'~')

    inner = recursive(obj.lhs,x,t,vars);
    rho = -inner;

elseif strcmp(obj.type,'<')

    lhs = recursive(obj.lhs,x,t,vars);
    rho = obj.rhs-lhs;

elseif strcmp(obj.type,'>')

    lhs = recursive(obj.lhs,x,t,vars);
    rho = lhs - obj.rhs;

elseif strcmp(obj.type,'<=')

    lhs = recursive(obj.lhs,x,t,vars);
    rho = obj.rhs-lhs + eps;

elseif strcmp(obj.type,'>=')

    lhs = recursive(obj.lhs,x,t,vars);
    rho = lhs - obj.rhs + eps;

elseif strcmp(obj.type,'+')

    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    rho = lhs + rhs;

elseif strcmp(obj.type,'-')

    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    rho = lhs - rhs;

elseif strcmp(obj.type,'*')

    lhs = obj.lhs;
    rhs = recursive(obj.rhs,x,t,vars);

    rho = lhs * rhs;

elseif strcmp(obj.type,'true')

    rho = inf*ones(size(t));

elseif strcmp(obj.type,'false')

    rho = -inf*ones(size(t));

elseif strcmp(obj.type,'variable')
    var_idx = strcmp(vars,obj.var);
    rho = x(:,var_idx);

elseif strcmp(obj.type,'finally')

    lhs = recursive(obj.lhs,x,t,vars);
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
    lhs = recursive(obj.lhs,x,t,vars);
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

    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    ind = find(t >= obj.from & t <= obj.to);
    rho=NaN(size(t)); cnt = 1;

    %check that time from is less than final time for traj
    if t(end) >= obj.from
        inner = NaN(size(rhs));
        for k = 1:length(inner)
            inner(k) = max(rhs(k),max(lhs(1:k)));
        end
        while ~isempty(ind) && ind(1) <= length(t)
            ind = ind(ind <= length(t));
            rho(cnt) = min(inner(ind));
            cnt = cnt + 1; ind = ind + 1;
        end
    end

elseif strcmp(obj.type,'until')
    lhs = recursive(obj.lhs,x,t,vars);
    rhs = recursive(obj.rhs,x,t,vars);

    ind = find(t >= obj.from & t <= obj.to);
    rho=NaN(size(t)); cnt = 1;

    %check that time from is less than final time for traj
    if t(end) >= obj.from
        inner = NaN(size(rhs));
        for k = 1:length(inner)
            inner(k) = min(rhs(k),min(lhs(1:k)));
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

