function R = reachKoopman(A,B,g,kfModel)
% compute reachable set for Koopman linearized model

%setup
dt=kfModel.ak.dt;
R0=kfModel.R0;
U=kfModel.U;
tFinal=kfModel.T;
cpBool=kfModel.cpBool;

% compute initial set using Taylor model arithmetic
n = dim(R0);
tay = taylm(R0);
tay = g(tay);
R0 = polyZonotope(tay);
R0 = polyZonotope(R0.c,R0.G,[],R0.expMat(1:n,:));

% compute reachable set
t = 0:dt:ceil(tFinal/dt)*dt;

set = cell(length(t),1); time = set; zono = set;

set{1} = R0; time{1} = interval(-dt/2,dt/2);

for i = 1:length(t)-1
    % AH edit to check if system has external input
    if B
        if kfModel.pulseInput
            cp_U = U.*cpBool(i,:)';
        else
            cp_U=U;
        end
        set{i+1} = A*set{i} + B*zonotope(cp_U);
    else
        set{i+1} = A*set{i};
    end
    time{i+1} = time{i} + dt;
end

% compute ouput set
for i = 1:length(set)
    set{i} = project(set{i},1:n);
    zono{i} = zonotope(set{i});

end

R.set = set; R.time = time; R.zono = zono;
end
