function g = randomFourierFeatureObservables(numFeat,dim,l)
% generate Random Fourier Feature observables cos(w'*x + u)

    % generate random scales and offsets
    w = normrnd(0,l^2,numFeat,dim);
    u = 2*pi*rand(numFeat,1);
    
    % generate fourier transform observables
    g = @(x) sqrt(2)*cos(w*x + u);
end