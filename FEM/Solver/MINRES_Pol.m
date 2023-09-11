classdef MINRES_Pol < handle

    properties (Access = public)
        xPrevIt
        convergenceData
        currentIter
    end

    methods (Access = public)

        function obj = MINRES_Pol()
            obj.init();
        end


        function x_k = solve(obj,A,b)
            obj.currentIter = obj.currentIter + 1;
            if isempty(obj.xPrevIt)
                obj.xPrevIt = zeros(size(b,1),1);
            end

            n = size(A,1);
            k = 3000;
            maxiter = k;
            tol = 10^(-4);


            [Alpha, Beta, v_0, V, c, s, Gamma, Delta, Epsilon, phi, d, x0, t] = obj.initiateVariables(k, n, b, A);



            j = 1;

            x_prev = x0;

            if (norm(A*x_prev-b)>tol)

                [Alpha(j), Beta(j+1), V(:,j+1)] = obj.lanczosProc(A, V(:,j), v_0, Beta(j));

                deltaj_prev = Delta(j);

                Delta(j)        =  c*deltaj_prev + s*Alpha(j);
                Gamma(j)        =  s*deltaj_prev - c*Alpha(j);
                Epsilon(j+1)    =  s*Beta(j+1);
                Delta(j+1)      = -c*Beta(j+1);

                [c, s, Gamma(j)] = obj.computeH(Beta, Gamma(j), j);
                % R(:,j) is now completed
                phi_prev = phi;
                phi      = s*phi_prev;
                t(j)     = c*phi_prev;
                d(:,j)   = V(:,j)/Gamma(j);

                x1 = x_prev + t(j)*d(:,j);

                j = 2;
                x_prev = x1;

                [Alpha(j), Beta(j+1), V(:,j+1)] = obj.lanczosProc(A, V(:,j), V(:,j-1), Beta(j));

                deltaj_prev = Delta(j);

                Delta(j)        =  c*deltaj_prev + s*Alpha(j);
                Gamma(j)        =  s*deltaj_prev - c*Alpha(j);
                Epsilon(j+1)    =  s*Beta(j+1);
                Delta(j+1)      = -c*Beta(j+1);

                [c, s, Gamma(j)] = obj.computeH(Beta, Gamma(j), j);
                % R(:,j) is now completed
                phi_prev = phi;
                phi      = s*phi_prev;
                t(j)     = c*phi_prev;
                d(:,j)   = (V(:,j)-Delta(j)*d(:,j-1))/Gamma(j);

                x2 = x_prev + t(j)*d(:,j);


                xj = x2;


                j=3;
                res = abs(phi);
                %             for j = 3:k
                while (j<maxiter)&&(res>tol)
                    x_prev = xj;

                    [Alpha(j), Beta(j+1), V(:,j+1)] = obj.lanczosProc(A, V(:,j), V(:,j-1), Beta(j));

                    deltaj_prev = Delta(j);

                    Delta(j)        =  c*deltaj_prev + s*Alpha(j);
                    Gamma(j)        =  s*deltaj_prev - c*Alpha(j);
                    Epsilon(j+1)    =  s*Beta(j+1);
                    Delta(j+1)      = -c*Beta(j+1);

                    [c, s, Gamma(j)] = obj.computeH(Beta, Gamma(j), j);            % R(:,j) is now completed

                    phi_prev    = phi;
                    phi         = s*phi_prev;
                    t(j)        = c*phi_prev;
                    d(:,j)      = (V(:,j)-Delta(j)*d(:,j-1)-Epsilon(j)*d(:,j-2))/Gamma(j);

                    xj = x_prev + t(j)*d(:,j);

                    res = abs(phi);
                    j=j+1;
                end
                x_prev = xj;
            end
            
            convIter = j-1
            x_k = x_prev;
            norm(A*x_k-b)/norm(b)
            obj.xPrevIt = x_k;
            obj.convergenceData(obj.currentIter) = convIter;
            if (obj.currentIter == 100)
                    stop = 1;
            end
        end
    end


    methods (Access = private)

        function init(obj)
            obj.currentIter = 0;
            obj.convergenceData = zeros(500,1);
        end

        function [alpha, beta_new, v_new] = lanczosProc(obj, A, vj, v_prev, betaj)
            z       = A*vj;
            alpha   = (z - v_prev*betaj)'*vj;
            w       = z - v_prev*betaj-vj*alpha;
            beta_new = norm(w);
            v_new    = w/beta_new;
        end

        function [Alpha, Beta, v_0, V, c, s, Gamma, Delta, Epsilon, phi, d, x0, t] = initiateVariables(obj, k, n, b, A)
            %             x0  = load('CounterPol\xNew.mat').xNew;
            x0 = zeros(size(b,1),1);
            % x0 = obj.xPrevIt;
            r0 = b-A*x0;

            Alpha   = zeros(k+1,1);
            Beta    = zeros(k+1,1);
            Beta(1) = norm(r0);
            V       = zeros(n,k);
            v_0     = zeros(n,1);


            V(:,1)  = r0/Beta(1);

            c       = -1;
            s       = 0;

            Gamma    = zeros(k,1);
            Delta    = zeros(k,1);
            Epsilon  = zeros(k,1);
            Gamma(1) = 0;

            phi = Beta(1);
            d   = zeros(n,k);


            t   = zeros(k,1);
        end


        function  [c, s, Gamma_j] = computeH(obj, Beta, Gamma_j, j)
            shi     = Gamma_j^2 + Beta(j+1)^2 + sign(Gamma_j)*Gamma_j*sqrt(Gamma_j^2+Beta(j+1)^2);
            c       = ((Beta(j+1)^2)/shi)-1;
            s       = -(Beta(j+1)*(Gamma_j + sign(Gamma_j)*sqrt(Gamma_j^2+Beta(j+1)^2)))/shi;

            Gamma_j = c*Gamma_j + s*Beta(j+1);
        end
    end
end