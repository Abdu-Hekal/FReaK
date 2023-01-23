function cost = compError(sys,g,X,U)
% compute the error between the Koopman linearized system and the
% measurements

    err = []; len = [];
    
    for r = 1:length(X)
        
       % simulation
       x = simulateKoopman(sys,g,X{r}(:,1),U{r});
       
       % error
       err = [err, sqrt(sum((X{r} - x).^2,1))];
       len = [len, sqrt(sum(X{r}.^2,1))];
    end
    
    % return different error measures
    cost.abs.max = max(err);
    cost.abs.mean = mean(err);
    cost.rel.max = max(err)/mean(len);
    cost.rel.mean = mean(err)/mean(len);
    cost.real.max = max(err./len);
    cost.real.mean = mean(err./len);
end