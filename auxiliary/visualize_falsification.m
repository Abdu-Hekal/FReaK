function visualize_falsification(crit_x, spec)
    %plot falsifying trace
    figure; hold on; box on;
    xlim([min(crit_x(:,1))-abs(0.2*min(crit_x(:,1))),max(crit_x(:,1))+abs(0.2*max(crit_x(:,1)))]); 
    ylim([min(crit_x(:,2))-abs(0.2*min(crit_x(:,2))),max(crit_x(:,2))+abs(0.2*min(crit_x(:,2)))])
    for i = 1:size(spec,1)
        if strcmp(spec(i,1).type,'unsafeSet')
            plot(spec(i,1).set, [1,2], 'FaceColor','red','FaceAlpha',.1)
        elseif strcmp(spec(i,1).type,'logic')
            phi = negationNormalForm(spec(i,1).set);
            plot_logic(phi);
        end
    end
    plot(crit_x(:,1),crit_x(:,2),'g');
end


% Auxiliary Functions -----------------------------------------------------

function plot_logic(phi)
    % convert logic equation to union of safe sets
    if ~phi.temporal
        eq = disjunctiveNormalForm(phi);
        safeSet = getClauses(eq,'dnf');
    
        for k = 1:length(safeSet)
            safeSet{k} = convert2set(safeSet{k});
        end
    
        % convert to a union of unsafe sets
        unsafeSet = safe2unsafe(safeSet);
        for k=1:length(unsafeSet)
            disp(unsafeSet{k}.P)
            plot(unsafeSet{k},[1,2], 'FaceColor','red','FaceAlpha',.1)
        end
    else
        plot_logic(phi.lhs);
        if phi.rhs
            plot_logic(phi.rhs);
        end
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