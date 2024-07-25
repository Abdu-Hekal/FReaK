function res = binarySplitInterval(interval, point)

dims=dim(interval);
splitInterval = split(interval,1:dims);
s1=splitInterval{1};
s2=splitInterval{2};

res=[];

for i=1:dims
    if contains(s1(i),point(i))
        res=vertcat(res,s1(i));
    elseif contains(s2(i),point(i))
        res=vertcat(res,s2(i));
    else
        soln.best.x(1,:)
        error('point must be within bounds of original interval')
    end
end
end
