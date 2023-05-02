classdef MINRES_RECYCLING < handle

    
    properties (Access = public)
        xPrevIt
        u_old
    end
    
    methods (Access = public)
        function obj = MINRES_RECYCLING()
            obj.init();
        end
        
        function x_k = solve(obj,A,b)

            n     = size(A,1);
            e1    = zeros(n,1);
            e1(1) = 1;
            k     = 5000; %max iterations
            s     = 100; %recycled vectors
            tol   = 10^(-4);

            [Alpha, Beta, v_0, V, c, s, Gamma, Delta, Epsilon, phi, d, x0, t] = obj.initiateVariables(k, n, b, A);

                x0 = obj.xPrevIt;
                r0 = b-A*x0;
                j = 1;


                %(altres, minres init)

                
            if exist(Yk,'var') == 1 %check if Yk is defined
                [Q,R] = qr(A*Yk,0);
                Ck = Q;
                Uk = Yk*inv(R);

                xRecycling = xPrev + Uk*Ck'*r0;
                rRecycling = r0 - Ck*Ck'*r0;


            else %Yk is not defined
                
                Beta(1) = norm(r0);

                V(:,1)  = r0/Beta(1);
                c0      = Beta(1)*e1;

                [rNew, xNew, Hk, Vk+1] = MINRES();

                Pk = createPk(n,s);
                Yk = V*Pk;
                [Q,R] = qr(Hk*Pk,0);
                Ck = Q;
                Uk = Yk*inv(R);
            end

            V(:,1)  = r0/norm;
            Anew = (I-Ck*Ck')*A;           
            [H, V] = obj.first2lanczos(Anew, V(:,j), v_0, Beta(j));
            j = 3;

            while (j<maxiter)&&(res>tol)
                [Alpha(j), Beta(j+1), V(:,j+1)] = lanczosProc(A, V(:,j), V(:,j-1), Beta(j));
                
                

                z = A*V;
                [rNew, xNew H, V, B] = MINRES();%k-s iteracions










                Pk = createPk(n,s);
                Yk = V*Pk;

                j = j+1;
            end









        end
    end

    methods (Access = private)
        function init(obj)
            obj.xPrevIt = zeros(20200,1);            
        end

        function [H, V] = obj.first2lanczos(obj, A, V(:,j), v_0, Beta(j))
            j = 1;
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

            H = zeros(k,k);
            H(1,1) = Alpha(1);
            H(2,2) = Alpha(2);
            H(2,1) = Beta(2);
            H(1,2) = Beta(2);
            H(2,3) = Beta(3);
            H(3,2) = Beta(3);
        end

        function [alpha, beta_new, v_new] = lanczosProc(obj, A, vj, v_prev, betaj)
            z       = A*vj;
            alpha   = (z - v_prev*betaj)'*vj;
            w       = z - v_prev*betaj-vj*alpha;
            beta_new = norm(w);
            v_new    = w/beta_new;
        end

        function Pk = createPk(obj,n,s)
                Pk = zeros(n,s);                
                Pk(1:s,1:s) = eye(s);
                Pk = Pk(randperm(n),:);
        end
        function [Alpha, Beta, v_0, V, c, s, Gamma, Delta, Epsilon, phi, d, x0, t] = initiateVariables(obj, k, n, b, A)
            %             x0  = load('CounterPol\xNew.mat').xNew;
            %             x0 = zeros(20200,1);
            x0 = obj.xPrevIt;
            r0 = b-A*x0;

            Alpha   = zeros(k+1,1);
            Beta    = zeros(k+1,1);
            Beta(1) = norm(r0);
            V       = zeros(n,k);
            v_0     = zeros(n,1);


%             V(:,1)  = r0/Beta(1);

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
    end
end

