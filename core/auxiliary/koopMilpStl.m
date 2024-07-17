function [F,P,O] = koopMilpStl(phi,kList,kMax,timePoints,var,M,normz,hardcoded,offsetMap)
% koopMilpStl - constructs constraints in YALMIP that compute
%             the robustness of satisfaction for stl specification phi
%
% Syntax:
%   [F,P,O] = koopMilpStl(phi,kList,kMax,timePoints,var,M,normz,hardcoded,offsetMap)
%
% Input:
%       phi:    an STLformula
%       kList:  a list of time steps at which the formula is to be enforced
%       kMAx:   the length of the trajectory
%       timePoints:  time points for analyzing robustness
%       var:    a dictionary mapping strings to variables
%       M:   	a large positive constant used for big-M constraints
%       normz:  normalization values based on boundaries of reachable sets
%       hardcoded: bool set to false if using optimizer object otherwise
%       true.
%       offsetMap: map that defines predicates to offset, key:number of
%       predicate, value:offset value
%
% Output:
%       F:  YALMIP constraints
%       P:  a struct containing YALMIP decision variables representing
%           the quantitative satisfaction of phi over each time step in
%           kList
%       O:  a struct containing YALMIP parameter variables representing the
%       offset for each inequality.
%
% :copyright: TBD
% :license: TBD

if (nargin==4)
    M = 1000;
end

F = [];
P = [];
O = [];

a = phi.from;
b = phi.to;

if isempty(a)
    a=0;
else
    % find idx of time point closest to a (rounded down, i.e. overapprox time horizon)
    a=find(timePoints <= a, 1, 'last')-1;
    %if all timePoints larger than a, then a is first timePoint
    if isempty(a)
        a=0;
    end
end
if isempty(b)
    b = kMax;
else
    % find idx of time point closest to b (rounded up, i.e. overapprox time horizon)
    b=find(timePoints >= b, 1, 'first')-1;
    %if all timePoints smaller than b, then b is last timePoint
    if isempty(b)
        b=numel(timePoints)-1;
    end
end

switch (phi.type)

    case {'<', '>', '<=', '>='} %predicate
        [F,P,O] = pred(phi,kList,var,normz,hardcoded,offsetMap);

    case '~'
        [Frest,Prest,Orest] = koopMilpStl(phi.lhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fnot, Pnot] = not(Prest);
        F = [F, Fnot, Frest];
        P = Pnot;
        O = [O, Orest];

    case '|'
        [Fdis1,Pdis1,O1] = koopMilpStl(phi.lhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fdis2,Pdis2,O2] = koopMilpStl(phi.rhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [For, Por] = or([Pdis1;Pdis2],M);
        F = [F, For, Fdis1, Fdis2];
        P = Por;
        O = [O, O1, O2];

    case '&'
        [Fcon1,Pcon1,O1] = koopMilpStl(phi.lhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fcon2,Pcon2,O2] = koopMilpStl(phi.rhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fand, Pand] = and([Pcon1;Pcon2],M);
        F = [F, Fand, Fcon1, Fcon2];
        P = Pand;
        O = [O, O1, O2];

    case 'globally'
        kListAlw = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a): min(kMax,k + b)}, kList)));
        [Frest,Prest,Orest] = koopMilpStl(phi.lhs,kListAlw,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Falw, Palw] = always(Prest,a,b,kList,kMax,M);
        F = [F, Falw];
        P = [Palw, P];
        F = [F, Frest];
        O = [O, Orest];

    case 'finally'
        kListEv = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a): min(kMax,k + b)}, kList)));
        [Frest,Prest,Orest] = koopMilpStl(phi.lhs,kListEv,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fev, Pev] = eventually(Prest,a,b,kList,kMax,M);
        F = [F, Fev];
        P = [Pev, P];
        F = [F, Frest];
        O = [O, Orest];

    case 'until'
        [Fp,Pp,O1] = koopMilpStl(phi.lhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Fq,Pq,O2] = koopMilpStl(phi.rhs,kList,kMax,timePoints, var,M,normz,hardcoded,offsetMap);
        [Funtil, Puntil] = until(Pp,Pq,a,b,kList,kMax,M);
        F = [F, Funtil, Fp, Fq];
        P = Puntil;
        O = [O, O1, O2];
end
end

function [F,z,varOffset] = pred(phi,kList,var,normz,hardcoded,offsetMap)
% Enforce constraints based on predicates
% var is the variable dictionary
fnames = fieldnames(var);

for ifield= 1:numel(fnames)
    eval([ fnames{ifield} '= var.' fnames{ifield} ';']);
end

st=str(phi);
%replace ge and le with gt and lt
st = replace(st,'<=','<');
st = replace(st,'>=','>');
if hardcoded
    global vkmrCount %globl count to track wihch subpred to offset in milp
    vkmrCount=vkmrCount+1; %increase count as pred found
    if isConfigured(offsetMap) && isKey(offsetMap,vkmrCount)
        robOffset=offsetMap(vkmrCount);
    else
        robOffset=0;
    end
    if contains(st, '<')
        tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
        st = ['-(' tokens{1}{1} '- (' tokens{1}{2} ')+' num2str(robOffset) ')'];
    end
    if contains(st, '>')
        tokens = regexp(st, '(.+)\s*>\s*(.+)','tokens');
        st= [ '(' tokens{1}{1} ')-(' tokens{1}{2} ')+' num2str(robOffset)];
    end
    varOffset = {}; %offset is hardcoded and so no spdvar.
else %use optimizer object
    varOffset=sdpvar(1,1); %setup rob offset
    if contains(st, '<')
        tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
        st = ['-(' tokens{1}{1} '- (' tokens{1}{2} ')+varOffset)'];
    end
    if contains(st, '>')
        tokens = regexp(st, '(.+)\s*>\s*(.+)','tokens');
        st= [ '(' tokens{1}{1} ')-(' tokens{1}{2} ')+varOffset' ];
    end
end

t_st = regexprep(st,'x(\w*)','x($1,t)'); %add time to state variables
t_st = regexprep(t_st,'u(\w*)','u($1,t)'); %add time to inputs
t_st = regexprep(t_st,'\<t\>',sprintf('%d:%d', kList(1),kList(end))); %add time points

%AH: Normalize
if isempty(normz)
else
    matches = regexp(t_st, 'x\((\d+),', 'tokens');
    tag = 'pred';
    normVal=0;
    for i=1:length(matches)
        idx = str2double(matches{i});
        tag = strcat(tag,',',num2str(idx));
        normVal = normVal + normz(idx);
    end
    t_st = strcat(t_st,'/',num2str(normVal));
end

try
    z_eval = eval(t_st);
end
zl = sdpvar(1,size(z_eval,2));
% F = [zl == z_eval]:tag;
F = zl == z_eval;
z = zl;
end


% BOOLEAN OPERATIONS

function [F,P] = and(p_list,M)
[F,P] = min_r(p_list,M);
end


function [F,P] = or(p_list,M)
[F,P] = max_r(p_list,M);
end


function [F,P] = not(p_list)
k = size(p_list,2);
P = sdpvar(1,k);
F = [P(:) == -p_list(:)];
end


% TEMPORAL OPERATIONS

function [F,P_alw] = always(P, a,b,kList,kMax,M)
F = [];
k = size(kList,2);
P_alw = sdpvar(1,k);
kListAlw = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a) : min(kMax,k + b)}, kList)));

for i = 1:k
    [ia, ib] = getIndices(kList(i),a,b,kMax);
    ia_real = find(kListAlw==ia);
    ib_real = find(kListAlw==ib);
    [F0,P0] = and(P(ia_real:ib_real)',M);
    F = [F;F0,P_alw(i)==P0];
end
end

%AH: fixed this function from blustl version
function [F,P_ev] = eventually(P, a,b,kList,kMax,M)
F = [];
k = size(kList,2);
P_ev = sdpvar(1,k);
kListEv = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a) : min(kMax,k + b)}, kList)));

for i = 1:k
    [ia, ib] = getIndices(kList(i),a,b,kMax);
    ia_real = find(kListEv==ia);
    ib_real = find(kListEv==ib);
    [F0,P0] = or(P(ia_real:ib_real)',M);
    F = [F;F0,P_ev(i)==P0];
end
end


function [F,P_until] = until(Pp,Pq,a,b,k,M)

F = [];
P_until = sdpvar(1,k);

for i = 1:k
    [ia, ib] = getIndices(i,a,b,kMax);
    F0 = [];
    P0 = [];
    for j = ia:ib
        [F1,P1] = until_mins(i,j,Pp,Pq,M);
        F0 = [F0, F1];
        P0 = [P0, P1];
    end
    [F4,P4] = max_r(P0);
    F = [F;F0,F4,P_until(i)==P4];
end

end


% UTILITY FUNCTIONS

function [F,P] = min_r(p_list,M)

k = size(p_list,2);
m = size(p_list,1);

P = sdpvar(1,k);
z = binvar(m,k);

F = [sum(z,1) == ones(1,k)];
for i=1:m
    F = [F, P <= p_list(i,:)];
    F = [F, p_list(i,:) - (1-z(i,:))*M <= P <= p_list(i,:) + (1-z(i,:))*M];
end
end

function [F,P] = max_r(p_list,M)

k = size(p_list,2);
m = size(p_list,1);

P = sdpvar(1,k);
z = binvar(m,k);

F = [sum(z,1) == ones(1,k)];
for i=1:m
    F = [F, P >= p_list(i,:)];
    %TODO: Our objective fcn is always to minimize robustness, thus binary
    %constraints here may not be needed (remove following line)
    F = [F, p_list(i,:) - (1-z(i,:))*M <= P <= p_list(i,:) + (1-z(i,:))*M];
end

end

function [F,P] = until_mins(i,j,Pp,Pq,M)
[F0,P0] = min_r(Pp(i:j)',M);
[F1,P] = min_r([Pq(j),P0],M);
F = [F0,F1];
end

function [ia, ib] = getIndices(i,a,b,k)
ia = min(k,i+a);
ib = min(k,i+b);
end



