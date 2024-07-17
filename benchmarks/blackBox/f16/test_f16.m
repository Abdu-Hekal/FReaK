init_f16;

altg = 4040;...
Vtg = 540;...  % Pass at Vtg = 540;    Fail at Vtg = 550;
phig = (pi/4);...       % Roll angle from wings level (rad)
thetag = -(pi/2)*0.8;...    % Pitch angle from nose level (rad)
psig = -(pi/4);             % Yaw angle from North (rad)

T = 0.1;

[tout, yout] = orig_run_f16(altg, Vtg, phig, thetag, psig, T);