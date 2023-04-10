function milp = KoopmanSetupReach(Sys)

Fstl=Sys.Fstl; Falpha=Sys.Falpha;
Alpha=Sys.Alpha; X=Sys.X; Pstl=Sys.Pstl;

%% Reachset constraints
Freach = [];
% Constraints for reachable set
for k=1:Sys.L+1
    % x = c + G * \alpha, 
    c = Sys.reachZonos{k}.center;
    G = Sys.reachZonos{k}.generators;

    Freach = [Freach, X(:,k) == (c+G*Alpha(1:size(G,2))')];
end


%% Objective function, minimize robustness of stl formula
obj = sum(sum(Pstl(:,1:end)));
milp = optimizer([Fstl, Falpha, Freach],obj,Sys.solver_options,[], {Alpha,X,Pstl});


