function [cost,g,A,B,C] = costFunPolynomial(x,X,U,dt,varargin)
% error for polynomial observables

    % parse input arguments
    f = []; outputMat = false;
    if nargin > 4 && ~isempty(varargin{1})
        f = varargin{1}; 
    end
    if nargin > 5 && ~isempty(varargin{2})
        outputMat = varargin{2};
    end

    % get properties
    numFeat = x(1);
    rank = x(2);
    rankOut = x(3);

    % generate observables
    n = size(X{1},1);
    g = polynomialObservables(numFeat,n,f);

    % identify system matrices
    sys = systemIdentification(X,U,g,dt,rank,outputMat,rankOut);

    % simulate Koopman linearized system and compute error
    A = sys.A; B = sys.B; C = sys.C;
    sys = linearSysDT(sys,dt);
    
    cost = 0;
    
    for r = 1:length(X)
        
       % simulation
       x = simulateKoopman(sys,g,X{r}(:,1),U{r});
       
       % error
       cost = cost + sum(sum((X{r} - x).^2,1));           
    end
end