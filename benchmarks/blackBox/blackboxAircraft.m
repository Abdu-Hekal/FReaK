function [ret]=blackboxAircraft(T,X,u)

B0=0.07351;
B1=-1.5E-3;
B2=6.1E-4;
C0=0.1667;
C1=0.109;
m=74E+3;
g=9.81;
S=158;
rho=0.3804;

% Interpolate input values at the current time T
input_idx = find(u(:, 1) <= T, 1, 'last'); % Find the index of the closest time point (piecewise constant)
if isempty(input_idx)
    input_idx = 1; % Use the first row if T is earlier than the first time point
end
u = u(input_idx, 2:end); %select corresponding input

mat1= [(-S * rho * B0 * X(1,1) * X(1,1))/(2 *m) - g * sin(pi*X(2,1)/180 );
       (S * rho * C0 * X(1,1)) / (2 * m) - g * cos(pi * X(2,1)/180)/X(1,1);
        X(1,1) * sin(pi * X(2,1)/180)];

mat2= [ u(1,1)/m; 0; 0];

mat3= [(-S * rho * X(1,1) * X(1,1))/(2*m) * (B1 * u(1,2) + B2 * u(1,2) * u(1,2));
       (S * rho * C1 )/(2*m)* X(1,1) * u(1,2);
       0];

ret = mat1 + mat2 + mat3;
    
 
    