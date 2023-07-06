function [offsetMap]=bReachCulprit3(critU,critX,t,xt,x0,A,B,g,R)
%function that gets idx of predicate (subformula) that is responsible for robustness value

offsetMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
idx=0;

u=critU(:,2:end)';
x = g(x0);
for i = 1:size(u,2)
    x = [x, A*x(:,end) + B*u(:,i)];
end
size(x)
koopTrace = [x,[u', zeros(size(u,2)-1,1)]];
interpCritX = interp1(xt,koopTrace,t,'pchip');
