function g = polynomialObservables(numFeat,n,varargin)
% generate polynomial observables

    % parse input arguments
    f = [];
    
    if nargin > 2 && ~isempty(varargin{1})
       f = varargin{1}; 
    end
    
    % generate initial variables
    x = sym('x',[n,1]); y = x;
    
    if ~isempty(f)
        y = f(x);
    end
    
    % generate observables
    cnt = 1; m = length(y);
    
    while true
        c = combinator(cnt,m);
        c = c(sum(c,2) <= cnt,:);
        cnt = cnt + 1;
        if size(c,1) > numFeat
           [~,ind] = sort(sum(c,2));
           c = c(ind(2:numFeat+1),:)-1; c(1:m,:) = eye(m); break;
        end
    end
    
    g = sym(zeros(numFeat,1));
    
    for i = 1:size(c,1)
        g(i) = prod(y.^(c(i,:)'));
    end
    
    g = matlabFunction(g,'Vars',{x});
end