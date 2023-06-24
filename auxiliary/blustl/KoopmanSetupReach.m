function Sys = KoopmanSetupReach(Sys)

alpha=Sys.alpha; x=Sys.x; L=Sys.L;
%% Reachset constraints
Freach = [];
% Constraints for reachable set
for k=1:L+1
    % x = c + G * \alpha, 
    c = Sys.reachZonos{k}.center;
    G = Sys.reachZonos{k}.generators;
    
    Freach = [Freach, x(:,k) == (c+G*alpha(1:size(G,2))')];
end

%% Reachset constraints
Sys.Freach=Freach;

