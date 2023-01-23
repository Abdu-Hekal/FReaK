function vol = computeVolume(sys,g,V,W,X,U,time)
% comptue the average volume of the reachable set

    vol = 0;

    % loop over all measurements
    for i = 1:length(X)
       
        % compute reachable set
        R = reachabilityAnalysis(sys,g,X{i}(:,1),U{i},time{i}(end),V,W);
        
        % compute volume
        vol_ = 0;
        
        for j = 1:length(R.timePoint.set)
           vol_ = vol_ + volume(interval(R.timePoint.set{j}));
        end
        
        vol = vol + vol_/length(R.timePoint.set);
    end

    vol = vol/length(X);
end