classdef koopmanRegressionLayer < nnet.layer.RegressionLayer
        
    properties
        A           % current system state matrix
        B           % current input matrix
        dt          % time step size
        n           % number of original system states
        m           % number of inputs
    end
 
    methods
        function layer = koopmanRegressionLayer(A,B,dt,n)           
            layer.A = A;
            layer.B = B;
            layer.dt = dt;
            layer.n = n;
            layer.m = size(B,2);
        end

        function loss = forwardLoss(layer, Y, T)
            X = T(1:layer.n,:);
            U = T(layer.n+1:layer.n+layer.m,:);
            Y = [X;Y];
            der = (Y(:,3:end) - Y(:,1:end-2))./(2*layer.dt);
            diff = der - [layer.A, layer.B]*[Y(:,2:end-1);U(:,2:end-1)];
            loss = sum(sum(abs(diff)));
        end
    end
end