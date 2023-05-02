classdef actualMINRES < handle

    
    properties (Access = public)
        xPrevIt
        Yk
    end
    




    methods (Access = public)
        function obj = actualMINRES()
            obj.init();
        end

        function x_k = solve(obj, A, b)

            n   = size(A,1);
            k   = 5000;
            s   = 100;
            tol = 10^(-4);

            e1    = zeros(n,1);
            e1(1) = 1;
            x0 = obj.xPrevIt;
            r0 = b - A*x0;

            if isempty(obj.Yk)

                beta1  = norm(r0);
                V(:,1) = r0/beta1;
                c0     = beta1*e1;                

                [xNew, rNew, V, H] = regularMINRES.solve(A, b);
                
                convIter = size(V,2);
                x_k = xNew;
                r_k = rNew;
                Pk = obj.createPk(convIter,s);
                obj.Yk = V*Pk;
                obj.xPrevIt = x_k;
            else
                tic
                qrA = A*obj.Yk;
                toc
                [Q, R] = qr(qrA, 0);
                Ck = Q;
                Uk = obj.Yk/R;
                x_start = obj.xPrevIt + Uk*Ck'*r0;
                r_start = r0-Ck*Ck'*r0;

                V(:,1) = r_start/norm(r_start);
                newA = (eye(20200)-Ck*Ck')*A;

                [xNew, rNew, V, H] = regularMINRES.solve(newA, b); %%triga molt, pk?
                x_k = xNew;
                B = Ck'*A*V;
                norms = sqrt(sum(Uk.^2,1));
                Dk = diag(1./norms);
                W = [Ck V];
                G = [Dk B; zeros(k) H];
                Pk = createPk(k,s);
                obj.Yk = V*Pk;

                [Q,R] = qr(G*Pk,0);
                Ck = W*Q;
                Uk = obj.Yk/R;
            end

        end







    end





    methods (Access = private)

        function init(obj)            
            obj.xPrevIt = zeros(20200,1);
        end

        function Pk = createPk(obj, k, s)
            Pk = zeros(k,s);
            Pk(1:s,1:s) = eye(s);
            Pk = Pk(randperm(k),:);
        end

    end




end

