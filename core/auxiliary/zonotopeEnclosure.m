function Z = zonotopeEnclosure(pZ,alpha)
% zonotopeEnclosure - enclose a polynomial zonotope by a zonotope
%
% Syntax:
%    Z = zonotopeEnclosure(pZ,alpha)
%
% Description:
%    Encloses a polynomial zonotope by a zonotope by linearizing the 
%    polynomial function around the linearization point for the factors 
%    alpha.
%
% Inputs:
%    pZ - polynomial zonotope
%    alpha - vector representing the linearization point for the factors
%
% Outputs:
%    Z - zonotope enclosing the polynomial zonotope
%
% Example:
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    Grest = [0;0.5];
%    expMat = [1 0 3;0 1 1];
%    pZ = polyZonotope(c,G,Grest,expMat);
% 
%    Z = zonotopeEnclosure(pZ,[0.5;0.5]);
%
%    plot(pZ);
%    hold on
%    plot(Z);


    % get linearization point
    f = pZ.c;

    for i = 1:size(pZ.G,2)
        f = f + pZ.G(:,i)*prod(alpha.^pZ.expMat(:,i));
    end

    % compute derivative for the polynomial function
    A = zeros(length(pZ.c),size(pZ.expMat,1));

    for i = 1:size(pZ.expMat,1)
        for j = 1:size(pZ.G,2)
            if pZ.expMat(i,j) > 0
                e = pZ.expMat(:,j);
                e(i) = e(i) - 1;
                A(:,i) = A(:,i) + (e(i)+1)*pZ.G(:,j)*prod(alpha.^e);
            end
        end
    end 

    f = f - A*alpha;
    pZlin = polyZonotope(f,A,[],eye(size(pZ.expMat,1)),pZ.id);

    % compute zonotope enclosure
    Ztemp = zonotope(exactPlus(pZlin,-1*pZ));

    Z = zonotope(f+center(Ztemp),[A,generators(Ztemp)]);
end