function [soln,falsified,robustness,Bdata,newBest]=checkFalsification(soln,x,u,t,specs,inputInterpolation,method,verb)
falsified=false; robustness=inf;
Bdata=NaN; newBest=false;
for ii=1:numel(specs)
    spec=specs(ii);
    % different types of specifications
    if strcmp(spec.type,'unsafeSet')
        falsified = any(spec.set.contains(x'));
    elseif strcmp(spec.type,'safeSet')
        falsified = ~all(spec.set.contains(x')); %check this
    elseif strcmp(spec.type,'logic')
        if ~isempty(u)
            interpU = interp1(u(:,1),u(:,2:end),t,inputInterpolation); %interpolate input at same time points as trajectory
        else
            interpU=u;
        end
        [Bdata,~,robustness] = bReachRob(spec,t,x,interpU');
        vprintf(verb,3,"robustness value: %.3f after %d simulations with method: %s \n",robustness,soln.sims,method)
        if robustness <= 0 %if less than or equal zero
            falsified=true;
        end
    end
    if robustness < soln.best.rob
        vprintf(verb,2,"new best robustness!: %.3f after %d simulations due to: %s \n",robustness,soln.sims,method)
        soln.best.rob=robustness;
        soln.best.x=x; soln.best.u=u; soln.best.t=t; soln.best.spec=spec;
        newBest=true; %found a new best soln
    end
    if falsified
        break;
    end
end
end