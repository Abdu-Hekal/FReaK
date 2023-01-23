function [W,V] = conformantSynthesisKoopman(sys,g,X,U,dt,mu,len)
% compute the set of disturbances as well as the set of measurement errors
% required to make the system conformant with the measurements

    % get output matrix and constant offset
    n = size(X{1},1); r = sys.dim;
    c = zeros(n,1); C = eye(n);
    
    if ~isempty(sys.c)
        c = sys.c;
    end
    if ~isempty(sys.C)
        C = sys.C;
    end

    % get state, input, and disturbance propagation matrices A, B, D
    if isa(sys,'linearSys')
        m = size(sys.B,2);
        sys = linearSysDT(linearSys(sys.A,[sys.B,eye(size(sys.A))],[],c),dt);
        A = sys.A; B = sys.B(:,1:m); D = sys.B(:,m+1:end); c = sys.c;
    else
        A = sys.A; B = sys.B; D = eye(size(A)); c = sys.c;
    end
    
    % initialization
    options = optimoptions('linprog','display','off');
    W = []; V = [];

    % loop over all measurements
    for i = 1:length(X)
       
        % split data into smaller chunks
        chunks = 1:len:size(X{i},2);
        if chunks(end) ~= size(X{i},2)
           chunks = [chunks size(X{i},2)];
        end
        
        x0 = g(X{i}(:,1)); R = x0;
        
        % loop over all data chunks
        for j = 1:length(chunks)-1
           
            X_ = X{i}(:,chunks(j)+1:chunks(j+1));
            U_ = U{i}(:,chunks(j):chunks(j+1)-1);
            
            % check if measurement is already contained in reachable set
            %     ==> no need to update uncertainty
            if i ~= 1 || j ~= 1
                
                points = zeros(n,len);
                
                if any(rad(W) > 0)
                
                    Wzono = zonotope(W); 

                    for k = 1:size(X_,2)
                        R = A*R + B*U_(:,k) + c + D*Wzono;
                        p = getClosestPoint(C*R,X_(:,k));
                        points(:,k) = X_(:,k) - p;
                    end

                    R = reduceUnderApprox(R,'sum',5);
                    
                else     
                    for k = 1:size(X_,2)
                        R = A*R + B*U_(:,k) + c + D*center(W);
                        points(:,k) = X_(:,k) - C*R;
                    end
                end
                
                if in(V,points)
                    if isa(R,'zonotope')
                        x0 = center(R);
                    else
                        x0 = R;
                    end
                    continue;
                end
            end
            
            % construct equality constraints H*x <= d for linear program
            H1 = []; H2 = []; d = []; x = x0; Dall = D;
            
            for k = 1:size(X_,2)
                H1 = [H1 zeros(size(H1,1),r);C*Dall];
                Dall = [A*Dall, D];
                H2 = blkdiag(H2,eye(n));
                x = A*x + B*U_(:,k) + c;
                d = [d; X_(:,k) - C*x];
            end
            
            H = [H1 -H1 H2 -H2];
            x0 = x;
            
            % construct objective function min_x f*x for linear program
            f = [mu*ones(2*size(H1,2),1); (1-mu)*ones(2*size(H2,2),1)];
            
            % solve linear program
            x = linprog(f,[],[],H,d,zeros(size(f)),[],options);
            
            % extract values for disturbances w_i and measurement err. v_i
            temp = blkdiag([eye(size(H1,2)) -eye(size(H1,2))], ...
                                     [eye(size(H2,2)) -eye(size(H2,2))])*x;
            w = temp(1:size(H1,2)); v = temp(size(H1,2)+1:end);
            w = reshape(w,[r,length(w)/r]);
            v = reshape(v,[n,length(v)/n]);
            
            % update sets for disturbances and measurement errors
            W = W | interval(min(w,[],2),max(w,[],2));
            V = V | interval(min(v,[],2),max(v,[],2));
            
            % initialize reachable set in the first iteration
            if i == 1 && j == 1
               R = x0;
            end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function p = getClosestPoint(R,p)
% estimate for the point closest to p in the zonotope R

    c = p - R.Z(:,1);
    alpha = sign(c'*R.Z(:,2:end));
    p = R.Z(:,1) + R.Z(:,2:end)*alpha';
end