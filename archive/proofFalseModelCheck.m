y = [linspace(0,1,11)',linspace(0,1,11)'];
t = linspace(0,1,11)';
x = stl('x',2);

% eq = globally(x(1) < 0.5,interval(0,1));
% res = checkStl(eq,y,t);
% assert(res==0)
% disp("checkStl correctly identifies that x(1) is not always greater than 0.1")
% res = modelCheckTrace(eq,y,t);
% assert(res==1)
% disp("modelCheckTrace fails")
% 
% eq = finally(x(1) > 1.1,interval(0,1));
% res = checkStl(eq,y,t);
% assert(res==0)
% disp("checkStl correctly identifies that x(1) is never greater than 1.1")
% res = modelCheckTrace(eq,y,t);
% assert(res==1)
% disp("modelCheckTrace fails")


%%conservativeness of checkstl
% eq = until(x(1) <= 0.5, x(2) >0.5, interval(0,1));
% res = checkStl(eq,y,t);
% assert(res==0)
% disp("checkStl is conservative, more suitable for reachability instead of falsification")

%%incorrect behaviour with nested functions
eq=finally(globally(x(2)>1.1,interval(0.1,0.3)),interval(0,1));
res = checkStl(eq,y,t);
spec=specification(eq,'logic');
[~,~,rob]=bReachRob(spec,t,y,[]);
assert(res==1)
disp('checkStl returns true for this case where it should be false, as x is never greater than 1')
assert(rob<0)
disp('breach correctly gives a negative value')

eq=globally(finally(x(2)>-0.1,interval(0.1,0.3)),interval(0,1));
res = checkStl(eq,y,t);
spec=specification(eq,'logic');
[~,~,rob]=bReachRob(spec,t,y,[]);
assert(res==0)
disp('checkStl returns false for this case where it should be true, as x is always greater than -0.1')
assert(rob>0)
disp('breach correctly gives a positive value')



