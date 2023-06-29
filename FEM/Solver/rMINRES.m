classdef rMINRES < handle

    properties        
        Uprev
        n
        maxIter
        nRecycling
        x0
        r0        
        A
        b 

        V
        B
        T
        S
        Phi
    end

    methods (Access = public)
        function x = solve(obj, A, b)
            obj.n           = size(b,1);
            obj.maxIter     = 5000;
            obj.nRecycling  = 25;
            tol             = 5e-5;

            obj.prepareProblem(A, b);
            j   = 0;
            x   = obj.x0;
            if (norm(obj.r0) > tol)
                [U, C, res, empty]          = obj.prepareRecycling(); 
                [G1,G2,g1,g2,d1,d2,z,beta]  = obj.prepareLoop(res, C);
                while ((j<obj.maxIter)&&(res>tol))||j<(obj.nRecycling+1)
                    j = j+1;
                    xOld = x;
                    [x, beta, G1, G2, d1, d2] = obj.computeNewX(j, C, beta, G1, G2, d1, d2, xOld);
                    [z, g1, g2]               = obj.computeNewZ(j, g1, g2, z);
                    res = norm(obj.Phi(j+1));
                end
                disp(j);
                x       = x + U*z;
                obj.x0  = x;                
                if empty == 1
                    obj.selectFirstRecycleSubspace(j);
                else
                    obj.selectRecycleSubspace(U,j,C);
                end         
            else
                x = obj.x0;
            end            
        end


        function obj = rMINRES()
            obj.init();
        end
    end

    
    methods (Access = private)
        function init(obj)
        end

        function prepareProblem(obj, A, b)
            if isempty(obj.x0)
                obj.x0 = zeros(obj.n, 1);
            end
            obj.r0 = b-A*obj.x0;
            obj.A = A;
            obj.b = b;
            obj.V   = zeros(obj.n, obj.maxIter+1);            
            obj.T   = zeros(obj.maxIter+1, obj.maxIter+1);
            obj.S   = zeros(obj.maxIter+1, obj.maxIter+1);
            obj.B   = zeros(obj.nRecycling, obj.maxIter+1);
            obj.Phi = zeros(obj.maxIter+1, 1);
        end

        function [U, C, res, empty] = prepareRecycling(obj)
            if isempty(obj.Uprev)
                empty       = 1;
                C           = zeros(obj.n,obj.nRecycling);
                U           = zeros(obj.n,obj.nRecycling);
                res         = norm(obj.r0);
                obj.V(:,1)  = obj.r0/res;
            else
                empty       = 0;
                t           = obj.A*obj.Uprev;
                [C, R]      = qr(t, 0);
                U           = obj.Uprev/R;
                rStart      = obj.r0-C*(C'*obj.r0);
                res         = norm(rStart);
                obj.V(:,1)  = rStart/res;
            end
        end

        function [G1,G2,g1,g2,d1,d2,z,beta] = prepareLoop(obj, res, C)
            beta   = norm(obj.r0);
            obj.Phi(1) = res;
            z = C'*obj.r0;
            G1 = [];
            G2 = [];
            d1 = 0;
            d2 = 0;
            g1 = 0;
            g2 = 0;
        end

        function [alphaNew, betaNew] = arnoldiStep(obj, j, C, vold, beta)
            v            = obj.V(:,j);
            a            = obj.A*v;  
            bj           = C'*a;
            w            = a - C*bj - vold*beta;
            alphaNew     = v'*w;
            w            = w - v*alphaNew;
            betaNew      = norm(w);
            obj.V(:,j+1) = w/betaNew;
            obj.B(:,j)   = bj;
        end
        
        function updateT(obj, j, alpha, beta)
            obj.T(j,j)     = alpha;
            obj.T(j,j+1)   = beta;
            obj.T(j+1,j)   = beta;
        end

        function [G1, G2] = updateS(obj, j, G1, G2)
            colS = obj.T(:,j);
            colS     = obj.firstUpdate(j, colS, G2);
            colS     = obj.secondUpdate(j, colS, G1);
            [c, s]   = obj.computeNewReflector(j, colS);
            colS     = obj.thirdUpdate(j, colS, c, s);
            obj.S(:,j) = colS;
            G2 = G1;
            G1 = [c s; s -c];
        end

        function updatedColumn = firstUpdate(obj, j, oldColumn, updaterOld2)
        t = oldColumn;
        if j > 2
            t(j-2:j-1) =  updaterOld2*t(j-2:j-1);
        end
        updatedColumn = t;
        end

        function updatedColumn = secondUpdate(obj, j, oldColumn, updaterOld1)
        t = oldColumn;
        if j > 1
            t(j-1:j) =  updaterOld1*t(j-1:j);
        end
        updatedColumn = t;
        end

        function updatedColumn = thirdUpdate(obj, j, oldColumn, c, s)
        t       = oldColumn;
        t(j)    = c*t(j)+s*t(j+1);
        t(j+1)  = 0;
        updatedColumn = t;
        end

        function [c, s] = computeNewReflector(obj, j, Scol)
            s1  = Scol(j);
            s2  = Scol(j+1);
            shi = s1^2 + s2^2 + sign(s1)*s1*sqrt(s1^2+s2^2);
            c   = ((s2^2)/shi)-1;
            s   = -(s2*(s1 + sign(s1)*sqrt(s1^2+s2^2)))/shi;
        end
 
        function updatePhi(obj, j, G1)
            PhiJ         = obj.Phi(j);
            obj.Phi(j)   = G1(1)*PhiJ;
            obj.Phi(j+1) = G1(2)*PhiJ;
        end

        function [d, d1, d2] = updateD(obj, j, d1, d2)
            colS = obj.S(:,j);
            v    = obj.V(:,j);
            sj   = colS(j);
            if j == 1
                s1 = 0;
                s2 = 0;
            elseif j == 2
                s1 = colS(j-1);
                s2 = 0;
            else
                s1 = colS(j-1);
                s2 = colS(j-2);
            end
            d  = (1/sj)*(v-d1*s1-d2*s2);
            d2 = d1;
            d1 = d;
        end




        function normalizer = buildNormalizerMatrix(obj, matrix)
            cols = size(matrix,2);
            normalizer = zeros(cols);
            for i = 1:cols
                normalizer(i,i) = 1/norm(matrix(:,i));
            end
        end

        function selectFirstRecycleSubspace(obj,j)
            m       = obj.nRecycling;
            Tbar    = obj.T(1:j+1,1:j);
            Tsquare = obj.T(1:j,1:j);
            Tend    = Tbar(j+1,j);
            vHat    = obj.V(:,1:j);
            ej      = zeros(j,1);
            ej(1)   = 1;
            LHS = (Tsquare+((Tend^2)*(Tsquare*(ej*ej'))));
            [eigVectors, eigValues]  = eig(LHS);
            eigValues       = diag(eigValues);
            [~, indices]    = mink(eigValues, m);
            obj.Uprev = vHat*eigVectors(:,indices);
        end

        function selectRecycleSubspace(obj,U,j,C)
            m              = obj.nRecycling;
            N              = obj.buildNormalizerMatrix(U);
            uTilde         = U*N;
            vHat           = [uTilde obj.V(:,1:j)];
            [LHS,RHS]       = obj.prepareGenEigenProblem(uTilde,N,j,C,m);
            [eigVectors,~]  = eigs(LHS, RHS, m, 'smallestabs','Tolerance', 1e-3);
            obj.Uprev = vHat*eigVectors;
        end

        function [LHS, RHS] = prepareGenEigenProblem(obj,uTilde,N,j,C,m)
            obj.B   = obj.B(:,1:j);
            Tfull       = obj.T(1:j+1,1:j);            
            NN = zeros(m);
            for i = 1:m
                NN(i,i) = N(i,i)^2;
            end
            NB = zeros(m,j);
            for i = 1:m
                NB(i,:) = obj.B(i,:)*N(i,i);
            end
            BN = zeros(j,m);
            for i = 1:m
                BN(:,i) = obj.B(i,:)*N(i,i);
            end
            LHS = [NN NB;BN (obj.B'*obj.B+Tfull'*Tfull)];
            CU  = C'*uTilde;
            VU  = obj.V(:,1:j+1)'*uTilde;
            NCU = zeros(m);
            for i = 1:m
                NCU(i,:) = CU(i,:)*N(i,i);
            end
            RHS = [NCU zeros(m, j); (obj.B'*CU + Tfull'*VU) Tfull(1:j,:)'];
        end

        function [x, beta, G1, G2, d1, d2] = computeNewX(obj, j, C, beta, G1, G2, d1, d2, x0ld)
            if j == 1
                vOld = 0;
            else
                vOld = obj.V(:,j-1);
            end
            [alpha, beta]   = obj.arnoldiStep(j, C, vOld, beta);            
            obj.updateT(j, alpha, beta);
            [G1, G2]        = obj.updateS(j, G1, G2);
            obj.updatePhi(j, G1);
            [d, d1, d2]     = obj.updateD(j, d1, d2);
            x = x0ld + d*obj.Phi(j);
        end
        function [z, g1, g2] = computeNewZ(obj, j, g1, g2, z)
            Phij = obj.Phi(j);
            bj   = obj.B(:,j);
            colS = obj.S(:,j);
            sj   = colS(j);
            if j == 1
                s1 = 0;
                s2 = 0;
            elseif j == 2
                s1 = colS(j-1);
                s2 = 0;
            else
                s1 = colS(j-1);
                s2 = colS(j-2);
            end
            g  = (1/sj)*(bj-g1*s1-g2*s2);
            g2 = g1;
            g1 = g;
            incZ    = g*Phij;
            z       = z - incZ;
        end
    end
end