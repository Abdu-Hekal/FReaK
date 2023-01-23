init_set = [1.5,2.25];
T = 7;

[tout, yout] = run_vanderpol(init_set, T);
disp(yout)
plot(yout(:,1),yout(:,2),'g');