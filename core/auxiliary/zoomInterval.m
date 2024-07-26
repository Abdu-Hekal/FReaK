function res = zoomInterval(I, point, varargin)

assert(size(point,2)==1,'interval must be 1 dimensional along y')
assert(size(point,2)==1,'point must be 1 dimensional along y')
assert(size(I,2)==size(point,2),'interval and point must be same size')
assert(contains(I,point),'point must be within the interval')

if nargin > 2
    alpha=varargin{1};
else
    alpha=0.1;
end

r = (I.sup- I.inf)*alpha; %range of zoom on interval
res=interval(point-r/2,point+r/2); %zoom on interval
res = res & I; %force zoom interval to be within interval bounds

end
