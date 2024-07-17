function [tout, yout]=runAircraft(T,x0,u)

[tout, yout] = ode45(@(T, X) aircraftODE(T, X, u), [0, T], x0);

end


    
 
    