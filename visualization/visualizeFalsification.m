function visualizeFalsification(varargin)
% visualizeFalsification - visualize the falsifying trace (and possibly unsafe set)
%
% Syntax:
%    visualizeFalsification(critX, times, spec, plotVars)
%    visualizeFalsification(critX, times, spec, plotVars,xlabel)
%    visualizeFalsification(critX, times, spec, plotVars,xlabel,ylabel)
%
% Inputs:
%    critX - the falsifying trajectory
%    times - time points of the trajectory
%    spec - falsifying spec
%    plotVars - dimensions of variables to plot (matrix or scalar value).
%    If plotVars is scalar value, plot against time
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author:      Abdelrahman Hekal
% Written:      28-February-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
numArgs = length(varargin);
assert(numArgs>=4,'must have at least 4 args, critical trajectory, time points, falsifying spec and plot vars')
critX=varargin{1};
times=varargin{2};
spec=varargin{3};
plotVars=varargin{4};

for i = 1:size(spec,1)
    figure; hold on; box on;

    if ~any(size(plotVars)>[1,1]) %singular plot var, plot against time
        plot(times,critX(:,plotVars),'b','DisplayName','falsifying traj');
    else
        xlim([min(critX(:,1))-abs(0.2*min(critX(:,1))),max(critX(:,1))+abs(0.2*max(critX(:,1)))]);
        ylim([min(critX(:,2))-abs(0.2*min(critX(:,2))),max(critX(:,2))+abs(0.2*min(critX(:,2)))]);
        if strcmp(spec(i,1).type,'unsafeSet')
            plot(spec(i,1).set, [1,2], 'FaceColor','red','FaceAlpha',.1,'DisplayName','unsafe set')
            %plot falsifying trace
            plot(critX(:,plotVars(1)),critX(:,plotVars(2)),'b','DisplayName','falsifying traj');
        elseif strcmp(spec(i,1).type,'logic')
            phi = spec(i,1).set;
            [unsafeCell, from, to] = phiToUnsafeSet(disjunctiveNormalForm(~phi),0,-1);
            from = find(times==from);
            if to==-1; to = size(times,1); else to=find(times==to); end
            plotUnsafeCell(unsafeCell, plotVars);


            if from ~= 1 %there exists trajectory before critical time period
                plot(critX(1:from,plotVars(1)),critX(1:from,plotVars(2)), 'Color' , 'g','LineStyle','--','LineWidth',1,'DisplayName','trajectory before');
            end
            plot(critX(from:to,plotVars(1)),critX(from:to,plotVars(2)),'Color' ,'g','LineWidth',1,'DisplayName','falsifying trajectory');
            if to ~= size(critX,1) %there exists trajectory after critical time period
                plot(critX(to:end,plotVars(1)),critX(to:end,plotVars(2)), 'Color' , 'g','LineStyle','--','LineWidth',1,'DisplayName','trajectory after');
            end
        end
    end
    l = legend;
    l.Location = 'northoutside';

    % Add axis titles
    if numArgs>4
        xlabel(varargin{5});
    end
    if numArgs>5
        ylabel(varargin{6});
    end
end
end


% -------------------------- Auxiliary Functions --------------------------

%FIXME: fix function to handle visualization of different stl formula
function [unsafeCell, from, to] = phiToUnsafeSet(phi, from, to)
if phi.from
    from = min(phi.from, from);
end
if phi.to
    to = max(phi.to, to);
end

unsafeCell = {};
% convert logic equation to union of safe sets
if ~phi.temporal
    unsafeSet = getClauses(phi,'dnf');
    for k = 1:length(unsafeSet)
        unsafeSet{k} = convert2set(unsafeSet{k});
    end

    % convert to a union of unsafe sets
    unsafeCell= [unsafeCell,unsafeSet];
elseif strcmp(phi.type,'|')
    % OR implies that unsafe set is BOTH rhs & lhs
    [lhsUnsafeCell, from, to] = phiToUnsafeSet(phi.lhs, phi.from, phi.to);
    [rhsUnsafeCell, from, to] = phiToUnsafeSet(phi.rhs, phi.from, phi.to);
    unsafeCell= [unsafeCell, lhsUnsafeCell, rhsUnsafeCell];
elseif strcmp(phi.type,'&')
    % AND implies that unsafe set is intersection of rhs & lhs
    [lhsUnsafeCell,from, to] = phiToUnsafeSet(phi.lhs,phi.from, phi.to);
    [rhsUnsafeCell,from, to] = phiToUnsafeSet(phi.rhs,phi.from, phi.to);
    for kl=1:length(lhsUnsafeCell)
        for kr=1:length(rhsUnsafeCell)
            %convert to mptPolytope to compute intersections of halfpaces
            unsafeIntersect = mptPolytope(lhsUnsafeCell{kl}) & mptPolytope(rhsUnsafeCell{kr});
            % if isempty(unsafeIntersect) %requires optimization toolbox
            %     error('stl impossible to falsify, intersection of (negation) requirements with and is empty')
            % end
            unsafeCell{end+1} = unsafeIntersect;
        end
    end
elseif any(strcmp(phi.type,{'finally','globally'}))
    [lhsUnsafeCell,from, to] = phiToUnsafeSet(phi.lhs,phi.from, phi.to);
    unsafeCell= [unsafeCell,lhsUnsafeCell];
elseif strcmp(phi.type,'next')
    [lhsUnsafeCell,from, to] = phiToUnsafeSet(phi.lhs,phi.from, phi.from);
    unsafeCell= [unsafeCell,lhsUnsafeCell];

    %FIXME: function throws error for plotting
    %     elseif strcmp(phi.type,'until')
    %         % AND implies that unsafe set is intersection of rhs & lhs
    %         [lhsUnsafeCell,from, to] = phiToUnsafeSet(phi.lhs,from, to);
    %         [rhsUnsafeCell,from, to] = phiToUnsafeSet(phi.rhs,from, to);
    %         for kl=1:length(lhsUnsafeCell)
    %             for kr=1:length(rhsUnsafeCell)
    %                 unsafeDiff = mldivide(lhsUnsafeCell{kl},rhsUnsafeCell{kr})
    %             end
    %             if ~isempty(unsafeDiff)
    %                 unsafeCell{end+1} = unsafeDiff;
    %             end
    %         end
else
    error("support for stl of type %s is not added\n",phi.type);
end
end

function list = safe2unsafe(sets)
% convert a safe set defined by the union of multiple polytopes to an
% equivalent union of unsafe sets

list = reverseHalfspaceConstraints(sets{1});

for i = 2:length(sets)

    tmp = reverseHalfspaceConstraints(sets{i});

    list_ = {};

    for j = 1:length(tmp)
        for k = 1:length(list)
            if isIntersecting(list{k},tmp{j})
                list_{end+1} = list{k} & tmp{j};
            end
        end
    end

    list = list_;
end
end

function res = reverseHalfspaceConstraints(poly)
% get a list of reversed halfspace constraints for a given polytope
res = {};
poly = mptPolytope(poly);

for i = 1:length(poly.P.b)
    res{end+1} = mptPolytope(-poly.P.A(i,:),-poly.P.b(i));
end
end

function plotUnsafeCell(unsafeCell,plotVars)
for k=1:length(unsafeCell)
    try
        plot(unsafeCell{k},plotVars, 'FaceColor','red','FaceAlpha',.05,'DisplayName','spec')
    catch
        disp("failed to plot spec")
    end
end
end