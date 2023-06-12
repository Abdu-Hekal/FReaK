function Sys = KoopmanSetupReach(Sys)

Alpha=Sys.Alpha; X=Sys.X; L=Sys.L;
%% Reachset constraints
Freach = [];
% Constraints for reachable set
for k=1:L+1
    % x = c + G * \alpha, 
    c = Sys.reachZonos{k}.center;
    G = Sys.reachZonos{k}.generators;

    Freach = [Freach, X(:,k) == (c+G*Alpha(1:size(G,2))')];
end

%% Reachset constraints
Sys.Freach=Freach;

