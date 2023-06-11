classdef rMINRES < handle

    properties        
        Uprev
        x0
        r0
        n
        convergenceData
        currentIter
    end

    methods (Access = public)
        function x = solve(obj, A, b)
            obj.currentIter = obj.currentIter + 1;
            obj.n       = size(b,1);
            maxIter     = 3000;
            nRecycling  = 400;
            tol         = 1e-3;
            obj.prepareProblem(A, b);
            V   = zeros(obj.n, maxIter);
            j   = 0;
            T   = zeros(maxIter, maxIter);
            S   = zeros(maxIter, maxIter);
            B   = zeros(nRecycling, maxIter);
            Phi     = zeros(maxIter, 1);
            x   = obj.x0;

            if (norm(obj.r0) > tol)
                if isempty(obj.Uprev)
                    res = norm(obj.r0);
                    Phi(1) = res;
                    V(:,1) = obj.r0/res;
                    beta   = res;
                    while ((j<maxIter)&&(res>tol))||j<(nRecycling+1)
                        j = j+1;
                        if j == 1
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, [], V(:,j), 0, beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), [], []);
                        elseif j == 2
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, [], V(:,j), V(:,j-1), beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), G1, []);
                        else
                            beta = Tpar(2);
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, [], V(:,j), V(:,j-1), beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), G1, G2);
                        end
                        modif    = Phi(j);
                        Phi(j)   = G1(1)*modif;
                        Phi(j+1) = G1(2)*modif;
                        if j == 1
                            d  = (1/S(j,j))*(V(:,j));
                            d1 = d;
                        elseif j == 2
                            d  = (1/S(j,j))*(V(:,j)-d1*S(j-1,j));
                            d2 = d1;
                            d1 = d;
                        else
                            d  = (1/S(j,j))*(V(:,j)-d1*S(j-1,j)-d2*S(j-2,j));
                            d2 = d1;
                            d1 = d;
                        end

                        incX    = d*Phi(j);
                        x       = x + incX;
                        res     = norm(Phi(j+1));
                    end
                    obj.x0 = x;


                    Tbar    = T(1:j+1,1:j);
                    T       = T(1:j,1:j);
                    Tend    = Tbar(j+1,j);
                    ej      = zeros(j,1);
                    ej(1)   = 1;
                    vHat    = V(:,1:j);

                    LHS = (T+((Tend^2)*(T*(ej*ej'))));                    
                    [eigVectors, eigValues]  = eig(LHS);
                    eigValues       = diag(eigValues);
                    [~, indices]    = mink(eigValues, nRecycling);
                    obj.Uprev = vHat*eigVectors(:,indices);
                else
                    [U, C] = obj.prepareRecycling(A);
                    rStart = obj.r0-C*(C'*obj.r0);
                    Phi(1)  = norm(rStart);
                    V(:,1) = rStart/Phi(1);
                    z = C'*obj.r0;
                    res = norm(obj.r0);
                    beta   = res;

                    while ((j<(maxIter-1))&&(res>tol))
                        j = j+1;
                        if j == 1
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, C, V(:,j), 0, beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), [], []);
                        elseif j == 2
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, C, V(:,j), V(:,j-1), beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), G1, []);
                        else
                            beta = Tpar(2);
                            [V(:,j+1), Tpar, beta] = obj.arnoldiStep(A, C, V(:,j), V(:,j-1), beta);
                            T = obj.updateT(j, Tpar, T);
                            [S(:,j), G1, G2]    = obj.computeScol(j, T(:,j), G1, G2);
                        end
                        modif    = Phi(j);
                        Phi(j)   = G1(1)*modif;
                        Phi(j+1) = G1(2)*modif;


                        bj = C'*(A*V(:,j));
                        B(:,j) = bj;
                        if j == 1
                            d  = (1/S(j,j))*(V(:,j));
                            d1 = d;
                            g  = (1/S(j,j))*(bj);
                            g1 = g;
                        elseif j == 2
                            d  = (1/S(j,j))*(V(:,j)-d1*S(j-1,j));
                            d2 = d1;
                            d1 = d;
                            g  = (1/S(j,j))*(bj-g1*S(j-1,j));
                            g2 = g1;
                            g1 = g;
                        else
                            d  = (1/S(j,j))*(V(:,j)-d1*S(j-1,j)-d2*S(j-2,j));
                            d2 = d1;
                            d1 = d;
                            g  = (1/S(j,j))*(bj-g1*S(j-1,j)-g2*S(j-2,j));
                            g2 = g1;
                            g1 = g;
                        end
                        incX    = d*Phi(j);
                        incZ    = g*Phi(j);
                        x       = x + incX;
                        z       = z - incZ;
                        res     = norm(Phi(j+1));
                    end
                    disp(j);
                    x = x + U*z;
                    obj.x0 = x;


                    N       = obj.buildNormalizerMatrix(U);
                    uTilde  = U*N;
                    B       = B(:,1:j);
                    T       = T(1:j+1,1:j);
                    vHat    = [uTilde V(:,1:j)];
                    wHat    = [C V(:,1:j+1)];
                    
                    NN = zeros(nRecycling);
                    for i = 1:nRecycling
                        NN(i,i) = N(i,i)^2;
                    end
                    NB = zeros(nRecycling,j);
                    for i = 1:nRecycling
                        NB(i,:) = B(i,:)*N(i,i);
                    end
                    BN = zeros(j,nRecycling);
                    for i = 1:nRecycling
                        BN(:,i) = B(i,:)*N(i,i);
                    end
                    LHS = [NN NB;BN (B'*B+T'*T)];

                    CU = C'*uTilde;
                    VU = V(:,1:j+1)'*uTilde;
                    NCU = zeros(nRecycling);
                    for i = 1:nRecycling
                        NCU(i,:) = CU(i,:)*N(i,i);
                    end
                    RHS = [NCU zeros(nRecycling, j); (B'*CU + T'*VU) T(1:j,:)'];




                    % Gm      = [N B; zeros(j+1, nRecycling) T];
                    % LHS     = Gm'*Gm;
                    % RHS     = Gm'*(wHat'*vHat);

                    % [eigVectors, eigValues]  = eig(LHS,RHS);      
                    % eigValues       = diag(eigValues);
                    % [~, indices]    = mink(eigValues, nRecycling);
                    % eigVectors = eigVectors(:,indices);
                    
                    % [eigVectors,~] = eigs(LHS,RHS,nRecycling,'smallestabs');
    

                    [eigVectors,~] = eigs(LHS, RHS, nRecycling, 'smallestabs','Tolerance', 1e-3);


                    obj.Uprev = vHat*eigVectors;


                    % obj.Uprev       = vHat(:,indices);


                    % [~, eigValues]  = eigs(LHS,RHS,nRecycling,'smallestreal');
                    % obj.Uprev = allVectors(:,randperm(size(allVectors,2),nRecycling));%RANDOM

                    


                end
            else
                x = obj.x0;
            end
            obj.convergenceData(obj.currentIter) = j;
            if (obj.currentIter == 100)
                stop = 1;
            end
        end

        function obj = rMINRES()
            obj.init();
        end
    end

    
    methods (Access = private)
        function init(obj)
            obj.currentIter = 0;
            obj.convergenceData = zeros(500,1);
        end

        function prepareProblem(obj, A, b)
            if isempty(obj.x0)
                obj.x0 = zeros(obj.n, 1);
            end
            obj.r0 = b-A*obj.x0;           
        
        end

        function [U, C] = prepareRecycling(obj, A)
            t       = A*obj.Uprev; 
            [C, R]  = qr(t, 0);
            U       = obj.Uprev/R;
        end

        function [vNew, Tpar, betaNew] = arnoldiStep(obj, A, C, v, vold, beta)
            a           = A*v;   
            if isempty(C)
                C = 0;
                b = 0;
            else
                b = C'*a;
            end
            w           = a - C*b - vold*beta;
            alphaNew    = v'*w;
            w           = w - v*alphaNew;
            betaNew     = norm(w);
            vNew        = w/betaNew;
            Tpar    = [0, 0];
            Tpar(1) = alphaNew;
            Tpar(2) = betaNew;
        end

        function T = updateT(obj, j, Tpar, T)
            T(j,j)      = Tpar(1);
            T(j+1,j)    = Tpar(2);
            T(j,j+1)    = Tpar(2);
        end
        function [Scol, G1, G2] = computeScol(obj, j, Tcol, G1, G2)
            Scol = zeros(size(Tcol,1),1);
            if isempty(G1)
                Scol = Tcol;
            elseif isempty(G2)
                t = Tcol;
                c = G1(1);
                s = G1(2);                
                Scol(j-1) = c*t(j-1) + s*t(j);
                Scol(j)   = s*t(j-1) - c*t(j);
                Scol(j+1) = t(j+1);
            else
            c = G2(1);
            s = G2(2);
            Scol(j-2)   = s*Tcol(j-1);
            Scol(j-1)   = -c*Tcol(j-1);
            Scol(j:j+1) = Tcol(j:j+1);
            t = Scol;

            c = G1(1);
            s = G1(2);
            Scol(j-2) = t(j-2);
            Scol(j-1) = c*t(j-1) + s*t(j);
            Scol(j)   = s*t(j-1) - c*t(j);
            end

            [c, s]  = obj.computeGivensRotation(Scol(j:j+1));
            Scol(j)   = c*Scol(j) + s*Scol(j+1);
            Scol(j+1) = 0;
            if isempty(G1)
                G2 = [];
                G1 = [c s];                
            else
                G2 = G1;
                G1 = [c s];
            end
        end

        function [c, s] = computeGivensRotation(obj, Scol)
            s1 = Scol(end-1);
            s2 = Scol(end);

            shi = s1^2 + s2^2 + sign(s1)*s1*sqrt(s1^2+s2^2);
            c   = ((s2^2)/shi)-1;
            s   = -(s2*(s1 + sign(s1)*sqrt(s1^2+s2^2)))/shi;
        end

        function normalizer = buildNormalizerMatrix(obj, matrix)
            cols = size(matrix,2);
            normalizer = zeros(cols);
            for i = 1:cols
                normalizer(i,i) = 1/norm(matrix(:,i));
            end
        end

    end
end