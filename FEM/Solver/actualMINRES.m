classdef actualMINRES < handle

    
    properties (Access = public)
        xPrevIt
        rPrevIt
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
                v1 = r0/beta1;
                c0     = beta1*e1;                

                [xNew, rNew, V, H, convIter] = regularMINRES.solve(A, b, obj.xPrevIt, v1, r0);
                V_noLastUpdate = V(:,(1:convIter));
                x_k = xNew;
                r_k = rNew;
                Vtilde = V_noLastUpdate;
                totalVectors = size(Vtilde,2);
                Pk = obj.createPk(totalVectors,s);
                obj.Yk = Vtilde*Pk;
                qrA = H*Pk;
                [Q,R] = qr(qrA,0);
                Uk = obj.Yk/R;
                obj.Yk = Uk;
                obj.xPrevIt = x_k;
            else
                tic
                qrA = A*obj.Yk;
                toc
                [Q, R] = qr(qrA, 0);
                Ck = Q;
                Uk = obj.Yk/R;
                x_known = Uk*(Ck'*r0);
                x_start = x0+x_known;
                r_start = r0-Ck*(Ck'*r0); %%equivalent a A*(x_start)-b

                v1 = r_start/norm(r_start);
                newA = (eye(20200)-Ck*Ck')*A;
                
                % newA = @(x) A*x-Ck*(Ck'*(A*x));
                
                [xNew, rNew, V, H, convIter] = regularMINRES.solve(newA, b, x_start, v1, r_start); %%triga molt, pk?
                
                V_noLastUpdate = V(:,(1:convIter));
                B = Ck'*A*V_noLastUpdate;

                norms = sqrt(sum(Uk.^2,1));
                Dk = diag(1./norms); %100x100
                Uktilde = Uk*Dk;
                Vtilde = [Uktilde V_noLastUpdate];
                W = [Ck V]; %20200x100 i 20200xconvIter+1
                G = [Dk B; zeros(convIter+1, s) H]; %100x100 / 100xconvIter // convIter+1x100 // convIter+1 x convIter 

                x_k = x_known + xNew;
                % r_k = r_start + rNew;              


                
                totalVectors = size(Vtilde,2);
                Pk = obj.createPk(totalVectors, s); %convIterx100
                
                obj.Yk = Vtilde*Pk; %20200x(convIter+s) x (convIter+s)x100
                
                qrA = G*Pk;
                [Q,R] = qr(qrA,0);
                Ck = W*Q;
                Uk = obj.Yk/R;
                obj.Yk = Uk;
                obj.xPrevIt = x_k;

            end

        end







    end





    methods (Access = private)

        function init(obj)            
            obj.xPrevIt = zeros(20200,1);
        end

        function Pk = createPk(obj, totalVectors, s)
            Pk = zeros(totalVectors,s);
            Pk(1:s,1:s) = eye(s);            
            Pk = Pk(randperm(totalVectors),:);
        end

    end




end

