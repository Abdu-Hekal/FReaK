function res = binarySplitInterval(I, point)

dims=dim(I);
splitInterval = split(I,1:dims);
s1=splitInterval{1};
s2=splitInterval{2};

res=[];

for i=1:dims
    if contains(s1(i),point(i))
        res=vertcat(res,s1(i));
    elseif contains(s2(i),point(i))
        res=vertcat(res,s2(i));
    else
        error('point must be within bounds of original interval')
    end
end
end
