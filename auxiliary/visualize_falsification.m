function visualize_falsification(crit_x, times, spec, plot_vars)
% visualize_falsification - visualize the falsifying trace (and possibly unsafe set)
%
% Syntax:  
%    visualize_falsification(crit_x, times, spec, plot_vars)
%
% Inputs:
%    crit_x - the falsifying trajectory 
%    times - time points of the trajectory
%    spec - falsifying spec 
%    plot_vars - dimensions of variables to plot (matrix or scalar value).
%    If plot_vars is scalar value, plot against time
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

    for i = 1:size(spec,1)
        figure; hold on; box on;
        if ~any(size(plot_vars)>[1,1]) %singular plot var, plot against time
            plot(times,crit_x(:,plot_vars),'b','DisplayName','falsifying traj');
        else
        xlim([min(crit_x(:,1))-abs(0.2*min(crit_x(:,1))),max(crit_x(:,1))+abs(0.2*max(crit_x(:,1)))]); 
        ylim([min(crit_x(:,2))-abs(0.2*min(crit_x(:,2))),max(crit_x(:,2))+abs(0.2*min(crit_x(:,2)))])
        if strcmp(spec(i,1).type,'unsafeSet')
            plot(spec(i,1).set, [1,2], 'FaceColor','red','FaceAlpha',.1,'DisplayName','unsafe set')
            %plot falsifying trace
            plot(crit_x(:,plot_vars(1)),crit_x(:,plot_vars(2)),'b','DisplayName','falsifying traj');
        elseif strcmp(spec(i,1).type,'logic')
            phi = negationNormalForm(spec(i,1).set);
            [unsafeCell, from, to] = phiToUnsafeSet(phi,0,-1);
            from = find(times==from);
            if to==-1; to = size(times,1); else to=find(times==to); end
            plotUnsafeCell(unsafeCell, plot_vars);
            plot(crit_x(1:from,plot_vars(1)),crit_x(1:from,plot_vars(2)), 'Color' , '#33F6FF','LineStyle','--','LineWidth',1.5,'DisplayName','trajectory before');
            plot(crit_x(from:to,plot_vars(1)),crit_x(from:to,plot_vars(2)),'b','LineWidth',1.5,'DisplayName','falsifying traj in spec range');
            plot(crit_x(to:end,plot_vars(1)),crit_x(to:end,plot_vars(2)), 'Color' , '#CA33FF','LineStyle','--','LineWidth',1.5,'DisplayName','trajectory after');
        end
        end
        l = legend;
        l.Location = 'best';
    end
end


% Auxiliary Functions -----------------------------------------------------

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
        eq = disjunctiveNormalForm(phi);
        safeSet = getClauses(eq,'dnf');
    
        for k = 1:length(safeSet)
            safeSet{k} = convert2set(safeSet{k});
        end
    
        % convert to a union of unsafe sets
        unsafeCell= [unsafeCell,safe2unsafe(safeSet)];
    elseif strcmp(phi.type,'&')
        % OR implies that unsafe set is BOTH rhs & lhs
        [lhsUnsafeCell, from, to] = phiToUnsafeSet(phi.lhs, phi.from, phi.to);
        [rhsUnsafeCell, from, to] = phiToUnsafeSet(phi.rhs, phi.from, phi.to);
        unsafeCell= [unsafeCell, lhsUnsafeCell, rhsUnsafeCell];
    elseif strcmp(phi.type,'|')
        % AND implies that unsafe set is intersection of rhs & lhs
        [lhsUnsafeCell,from, to] = phiToUnsafeSet(phi.lhs,phi.from, phi.to);
        [rhsUnsafeCell,from, to] = phiToUnsafeSet(phi.rhs,phi.from, phi.to);
        for kl=1:length(lhsUnsafeCell)
            for kr=1:length(rhsUnsafeCell)
                unsafeIntersect = lhsUnsafeCell{kl} & rhsUnsafeCell{kr};
                if ~isempty(unsafeIntersect)
                    unsafeCell{end+1} = unsafeIntersect;
                end
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
        fprintf("support for stl of type %s is not added\n",phi.type);
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

function plotUnsafeCell(unsafeCell,plot_vars)
    disp(length(unsafeCell))
    for k=1:length(unsafeCell)
        try
            plot(unsafeCell{k},plot_vars, 'FaceColor','red','FaceAlpha',.05,'DisplayName','spec')
        catch
            disp("failed to plot spec")
        end
    end
end