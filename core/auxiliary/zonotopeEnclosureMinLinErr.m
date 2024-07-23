function Z = zonotopeEnclosureMinLinErr(pZ,alpha,varargin)
% zonotopeEnclosureMinLinErr - enclose a polynomial zonotope by a zonotope
%
% Syntax:
%    Z = zonotopeEnclosureMinLinErr(pZ,alpha)
%    Z = zonotopeEnclosureMinLinErr(pZ,alpha,extreme)
%
% Description:
%    Encloses a polynomial zonotope by a zonotope by linearizing the 
%    polynomial function in such a way that the linearization error is 
%    minimized.
%
% Inputs:
%    pZ - polynomial zonotope
%    alpha - vector representing the linearization point for the factors
%    extreme - minimize lineariaztion error at extreme points only (extreme
%              = true) or also on other points (extreme = false = default)
%
% Outputs:
%    Z - zonotope enclosing the polynomial zonotope
%
% Example:
%    c = [1;-2;3];
%    G = [0 1 0;2 0 0;0 0 2];
%    Grest = [0;0.5;0.5];
%    expMat = [1 0 3;0 1 1];
%    pZ = polyZonotope(c,G,Grest,expMat);
%
%    Z = zonotopeEnclosureMinLinErr(pZ,[0.5;0.5]);
%
%    figure; hold on; plot(pZ,[1,2,3],'b'); plot(Z,[1,2,3],'r');

    % parse input arguments
    extreme = false;

    if nargin > 2
        extreme  = varargin{1};
    end

    % get linearization point
    f = pZ.c;

    for i = 1:size(pZ.G,2)
        f = f + pZ.G(:,i)*prod(alpha.^pZ.expMat(:,i));
    end

    % get factor values at which the linearization error is minimized
    n = size(pZ.expMat,1);
    I = interval(-ones(n,1),ones(n,1)); a = [];

    for i = 1:n
        tmp = zeros(n,1); tmp(i) = 1;
        a = [a,tmp,-tmp];
    end

    if 2^n < 100
        a = [a,vertices(I)];
    end

    if size(a,2) < 100
        if extreme
            a = [a,randPoint(I,100-size(a,2),'extreme')];
        else
            a = [a,randPoint(I,100-size(a,2))];
        end
    end

    % compute corresponding function values for the alpha values
    fval = zeros(length(pZ.c),n);

    for i = 1:size(a,2)
        fval(:,i) = pZ.c;
        fval(:,i) = fval(:,i) + sum(pZ.G.*prod(a(:,i).^pZ.expMat,1),2);
    end

    % linear regression to obtain a least-square linear fit 
    X1 = a - alpha;
    X2 = fval - f;

    [V,S,W] = svd(X1,'econ');

    m = min(size(S));
    V_ = V(:,1:m); S_ = S(1:m,1:m); W_ = W(:,1:m);

    A = X2*W_*diag(1./diag(S_))*V_';

    f = f - A*alpha;
    pZlin = polyZonotope(f,A,[],eye(size(pZ.expMat,1)),pZ.id);

    % compute zonotope enclosure
    Ztemp = zonotope(exactPlus(pZlin,-1*pZ));

    G = zeros(size(pZ.G));
    ind = find(sum(pZ.expMat,1) == 1);

     for i = 1:size(pZ.expMat,1)
        for j = ind
            if pZ.expMat(i,j) == 1
                G(:,j) = A(:,i);
            end
        end
    end

    Z = zonotope(f+center(Ztemp),[G,pZ.Grest]);
end