function sys = systemIdentification(X,U,g,dt,rank,varargin)
% identify the system matrices for the linear systems

    % parse input arguments
    outputMat = false; rankOut = []; offset = false;
 
    if nargin > 4 && ~isempty(varargin{1})
        outputMat = varargin{1};
    end
    
    if nargin > 5 && ~isempty(varargin{2})
        rankOut = varargin{2}; 
    end
    
    if nargin > 6 && ~isempty(varargin{3})
        offset = varargin{3}; 
    end

    % map measurements to the observables space
    X_ = cell(size(X)); numFeat = length(g(X{1}(:,1)));
    for r = 1:length(X)
       temp = X{r}; X_{r} = zeros(numFeat,size(temp,2));
       for h = 1:size(temp,2)
          X_{r}(:,h) = g(temp(:,h)); 
       end
    end

    % dynamic mode decomposition
    if offset

        [A,B,c] = dmd(X_,dt,U,'cont','svd',rank);
    else
        [A,B] = dmd(X_,dt,U,'cont','svd',rank);
        c = zeros(size(A,1),1);
    end
    
    % determine output matrix C
    n = size(X{1},1); numFeat = size(A,1);
    C = [eye(n), zeros(n,numFeat-n)];
    
    if outputMat
        sysDT = linearSysDT(linearSys(A,B,[],eye(numFeat)),dt);
        Y = [X{:}]; X_ = [];
        for i = 1:length(X)
           X_  = [X_, simulateKoopman(sysDT,g,X{i}(:,1),U{i})]; 
        end
        C = getOutputMatrix(Y,X_,rankOut);
    end

    % constrcut resulting linear system object
    sys = linearSys(A,B,c,C);
end