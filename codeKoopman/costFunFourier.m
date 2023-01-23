function [cost,g,A,B,C,c] = costFunFourier(x,X,U,dt,varargin)
% error for random fourier feature observables

    % parse input arguments
    outputMat = false; offset = false;
    
    if nargin > 4 && ~isempty(varargin{1})
        outputMat = varargin{1};
    end
    
    if nargin > 5 && ~isempty(varargin{2})
        offset = varargin{2}; 
    end

    % get properties
    numFeat = x(1);
    rank = x(2);
    l = x(3);
    rankOut = x(4);

    % generate observables
    n = size(X{1},1);
    
    if ~outputMat
        g_ = randomFourierFeatureObservables(numFeat-n,n,l);
        g = @(x) [x; g_(x)];
    else
        g = randomFourierFeatureObservables(numFeat,n,l);
    end


    % identify system matrices
    sys = systemIdentification(X,U,g,dt,rank,outputMat,rankOut,offset);

    % simulate Koopman linearized system and compute error
    A = sys.A; B = sys.B; C = sys.C; c = sys.c;
    sys = linearSysDT(sys,dt);
    
    err = []; len = [];
    
    for r = 1:length(X)
        
       % simulation
       x = simulateKoopman(sys,g,X{r}(:,1),U{r});
       
       % error
       err = [err, sqrt(sum((X{r} - x).^2,1))];
       len = [len, sqrt(sum(X{r}.^2,1))];
    end
    
    cost.abs.max = max(err);
    cost.abs.mean = mean(err);
    cost.rel.max = max(err)/mean(len);
    cost.rel.mean = mean(err)/mean(len);
    cost.real.max = max(err./len);
    cost.real.mean = mean(err./len);


    % return fourier feature observables as a symbolic function
    xSym = sym('x',[n,1]);
    g = g(xSym);

end