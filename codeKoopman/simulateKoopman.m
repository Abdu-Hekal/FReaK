function x = simulateKoopman(sys,g,x0,u)
% simulate a Koopman operator linearized system

    xtemp = g(x0);  
    x = zeros(size(sys.C,1),size(u,2)); 
    x(:,1) = sys.C*xtemp;

    for h = 1:size(u,2)
        xtemp = sys.A*xtemp + sys.B*u(:,h) + sys.c;
        x(:,h+1) = sys.C*xtemp;
    end
end