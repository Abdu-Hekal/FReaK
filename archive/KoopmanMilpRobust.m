function [F,P] = KoopmanMilpRobust(phi,kList,kMax,ts,var,M,normz,offsetMap)
% KoopmanMilpRobust  constructs MILP constraints in YALMIP that compute
%                  the robustness of satisfaction for specification phi
%
% Input:
%       phi:    an STLformula
%       kList:  a list of time steps at which the formula is to be enforced
%       kMAx:   the length of the trajectory
%       ts:     the interval (in seconds) used for discretizing time
%       var:    a dictionary mapping strings to variables
%       M:   	a large positive constant used for big-M constraints
%       normz:  normalization values based on boundaries of reachable sets
%       offsetMap: map that defines predicates to offset, key:number of
%       predicate, value:offset value
%
% Output:
%       F:  YALMIP constraints
%       P:  a struct containing YALMIP decision variables representing
%           the quantitative satisfaction of phi over each time step in
%           kList
%
% :copyright: TBD
% :license: TBD

if (nargin==4)
    M = 1000;
end

F = [];
P = [];

a = phi.from;
b = phi.to;

if isempty(a)
    a=0;
else
    a = max([0 floor(a/ts)]);
end
if isempty(b)
    b = kMax;
else
    b = ceil(b/ts);
end

switch (phi.type)

    case {'<', '>', '<=', '>='} %predicate
        [F,P] = pred(phi,kList,var,normz,offsetMap);

    case '~'
        [Frest,Prest] = KoopmanMilpRobust(phi.lhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Fnot, Pnot] = not(Prest);
        F = [F, Fnot, Frest];
        P = Pnot;

    case '|'
        [Fdis1,Pdis1] = KoopmanMilpRobust(phi.lhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Fdis2,Pdis2] = KoopmanMilpRobust(phi.rhs,kList,kMax,ts, var,M,normz,offsetMap);
        [For, Por] = or([Pdis1;Pdis2],M);
        F = [F, For, Fdis1, Fdis2];
        P = Por;

    case '&'
        [Fcon1,Pcon1] = KoopmanMilpRobust(phi.lhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Fcon2,Pcon2] = KoopmanMilpRobust(phi.rhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Fand, Pand] = and([Pcon1;Pcon2],M);
        F = [F, Fand, Fcon1, Fcon2];
        P = Pand;

    case 'globally'
        kListAlw = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a): min(kMax,k + b)}, kList)));
        [Frest,Prest] = KoopmanMilpRobust(phi.lhs,kListAlw,kMax,ts, var,M,normz,offsetMap);
        [Falw, Palw] = always(Prest,a,b,kList,kMax,M);
        F = [F, Falw];
        P = [Palw, P];
        F = [F, Frest];

    case 'finally'
        kListEv = unique(cell2mat(arrayfun(@(k) {min(kMax,k + a): min(kMax,k + b)}, kList)));
        [Frest,Prest] = KoopmanMilpRobust(phi.lhs,kListEv,kMax,ts, var,M,normz,offsetMap);
        [Fev, Pev] = eventually(Prest,a,b,kList,kMax,M);
        F = [F, Fev];
        P = [Pev, P];
        F = [F, Frest];

    case 'until'
        [Fp,Pp] = KoopmanMilpRobust(phi.lhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Fq,Pq] = KoopmanMilpRobust(phi.rhs,kList,kMax,ts, var,M,normz,offsetMap);
        [Funtil, Puntil] = until(Pp,Pq,a,b,kList,kMax,M);
        F = [F, Funtil, Fp, Fq];
        P = Puntil;
end
end

function [F,z] = pred(phi,kList,var,normz,offsetMap)
% Enforce constraints based on predicates
% var is the variable dictionary
fnames = fieldnames(var);

for ifield= 1:numel(fnames)
    eval([ fnames{ifield} '= var.' fnames{ifield} ';']);
end

global vkmrCount %globl count to track wihch subpred to offset in milp
vkmrCount=vkmrCount+1; %increase count as pred found
if isKey(offsetMap,vkmrCount)
    robOffset=offsetMap(vkmrCount);
else
    robOffset=0;
end

st=str(phi);
if contains(st, '<')
    tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
    st = ['-(' tokens{1}{1} '- (' tokens{1}{2} ')+' num2str(robOffset) ')'];
end
if contains(st, '>')
    tokens = regexp(st, '(.+)\s*>\s*(.+)','tokens');
    st= [ '(' tokens{1}{1} ')-(' tokens{1}{2} ')+' num2str(robOffset)];
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



