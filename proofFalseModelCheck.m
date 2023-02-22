y = linspace(0,1,100)';
t = linspace(0,1,100)';

x = stl('x',2); 
eq = globally(x(1) > 0.1,interval(0,1));
res = checkStl(eq,y,t);
assert(res==0)
disp("checkStl correctly identifies that x(1) is not always greater than 0.1")
res = modelCheckTrace(eq,y,t);
assert(res==1)
disp("modelCheckTrace fails")

eq = finally(x(1) > 1.1,interval(0,1));
res = checkStl(eq,y,t);
assert(res==0)
disp("checkStl correctly identifies that x(1) is never greater than 1.1")
res = modelCheckTrace(eq,y,t);
assert(res==1)
disp("modelCheckTrace fails")
