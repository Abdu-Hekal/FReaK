function maxTime = maxStlTime(stl)
maxTime=0;
timedOp = {'finally', 'globally', 'release', 'until'};
maxTime=recursive(stl,timedOp,maxTime);
end

function maxTime = recursive(stl,timedOp,maxTime)
if isa(stl,'stl')
    isMember = ismember(stl.type, timedOp);
    if isMember
        maxTime=maxTime+stl.to;
    end
    if ~isempty(stl.lhs)
        maxTime=recursive(stl.lhs,timedOp,maxTime);
    end
    if ~isempty(stl.rhs)
        maxTime=recursive(stl.rhs,timedOp,maxTime);
    end
end
end
