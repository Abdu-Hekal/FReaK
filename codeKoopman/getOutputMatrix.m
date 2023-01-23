function C = getOutputMatrix(Y,X,varargin)
% determine the best ouput matrix such that ||C*X - Y||_2 is minimized

    % parse input arguments
    rank = [];
    if nargin > 2 && ~isempty(varargin{1})
       rank = varargin{1}; 
    end

    % singular value decompostion
    [U,S,W] = svd(X,'econ');
    
    % reduce rank
    if ~isempty(rank)
        U = U(:,1:rank); S = S(1:rank,1:rank); W = W(:,1:rank);
    end
    
    % construct resulting C matrix
    C = Y*W*diag(1./diag(S))*U';
end