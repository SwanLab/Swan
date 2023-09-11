classdef actualMINRES < handle  

    properties (Access = public)
        x0
        r0
        Yk
        n
    end

    methods (Access = public)
        function obj = actualMINRES()
            obj.init();
        end

        function x_k = solve(obj, A, b)
            obj.n       = size(b,1);
            maxIter = 3000;
            omega   = 50;
            tol     = 10e-5;            
            obj.computeInitialValues(A, b);
            if isempty(obj.Yk)
                beta1  = norm(obj.r0);
                v1 = obj.r0/beta1;          
                [xNew, rNew, V, H, convIter] = regularMINRES.solve(A, b, obj.x0, v1, obj.r0, omega); %(!)
                V_noLastUpdate = V(:,(1:convIter));
                x_k = xNew;
                Vtilde = V_noLastUpdate;
                totalVectors = size(Vtilde,2);
                obj.Yk = Vtilde(:,randperm(totalVectors,omega));
                obj.x0 = x_k;
            else
                [Ck, R] = obj.computeQR(A);                
                r_start = obj.r0-Ck*(Ck'*obj.r0);
                res     = norm(r_start);

                [Alpha, Beta]   = obj.defineTVectors();
                [V, vTilde]     = obj.defineVVectors(maxIter);
                [bCols, bTilde] = obj.defineBVectors(maxIter, omega);
                S               = zeros(maxIter,maxIter);
                [bHat, yHat]    = obj.defineRecyclingVectors(Ck, R, res);

                
                V(:,1)  = r_start/res;
                x_old   = obj.x0;

                

                j = 1;

                vHat = A*V(:,j);
                vHat = vHat - Ck*(Ck'*vHat);
                bCols(:,j) = R\(Ck'*vHat); %b1

                Alpha(j)  = V(:,j)'*vHat; %ok
                vHat    = vHat - V(:,j)*Alpha(j); 
                Beta(j+1)   = norm(vHat); %ok
                V(:,j+1)      = vHat/Beta(j+1);

                %compute column 1 of T

                Tcol = [Alpha(j); Beta(j+1)];

                %apply givens rotation j-2, j-1 (only j-1 when j=0)

                S1   = Tcol;


                %compute G

                [c, s] = obj.computeGivensRotation(S1); %G1

                S(j,j) = c*S1(1)+s*S1(2);
                %S(j+1,j) = 0;


                modif       = yHat(j);
                yHat(j)     = c*modif;
                yHat(j+1)   = s*modif;

                vTilde(:,j) = (1/S(j,j))*(V(:,j)); %equivalent d
                bTilde(:,j) = (1/S(j,j))*(bCols(:,j));
                x_new      = x_old + vTilde(:,j)*yHat(j);
                bHat      = bHat - bTilde(:,j)*yHat(j);


                j = 2;
                x_old = x_new;

                vHat = A*V(:,j);
                vHat = vHat - Ck*(Ck'*vHat);
                bCols(:,j) = R\(Ck'*vHat);


                vHat        = vHat-V(:,j-1)*Beta(j);
                Alpha(j)    = V(:,j)'*vHat; 
                vHat        = vHat - V(:,j)*Alpha(j); 
                Beta(j+1)   = norm(vHat);
                V(:,j+1)    = vHat/Beta(j+1);

                %compute column 2 of T

                Tcol = [Beta(j); Alpha(j); Beta(j+1)];

                %apply givens rotation j-2, j-1

                % Gold_2 = Gold_1; %G0
                Gold_1 = [c s 0; s -c 0; 0 0 1]; %G1, same as current Gnew

                S2 = Gold_1*Tcol; %for j = 2, G0 doesnt affect                

                %compute G

                [c, s] = obj.computeGivensRotation(S2);    

                S(j-1,j) = S2(1);
                S(j,j)   = c*S2(2)+s*S2(3);
                %S(j+1,j) = 0;

                modif       = yHat(j);
                yHat(j)     = c*modif;
                yHat(j+1)   = s*modif;

                vTilde(:,j) = (1/S(j,j))*(V(:,j)-vTilde(:,j-1)*S(j-1,j)); %equivalent d
                bTilde(:,j) = (1/S(j,j))*(bCols(:,j)-bTilde(:,j-1)*S(j-1,j));
                x_new      = x_old + vTilde(:,j)*yHat(j);
                bHat      = bHat - bTilde(:,j)*yHat(j);

                j = 3;
                x_old = x_new;

                Gnew = [1 0 0; 0 c s; 0 s -c];

                while ((j<maxIter)&&(norm(yHat(j))>(tol*res)))||j<(omega+1)
                    vHat = A*V(:,j);
                    vHat = vHat - Ck*(Ck'*vHat);
                    bCols(:,j) = R\(Ck'*vHat);


                 [incX,incB,S,vTilde,bTilde,V,yHat,Alpha,Beta,Gold_1,Gnew] = obj.Lanczos(vHat,bCols,S,vTilde,bTilde,V,yHat,Alpha,Beta,Gold_1,Gnew,j);

                    x_new       = x_old + incX;
                    bHat        = bHat - incB;

                    x_old   = x_new;
                    j       = j+1;
                end
                x_known = obj.Yk*bHat;
                x_k     = x_new + x_known;
                norm(A*x_k-b)
                norm(A*x_new-b)
                norm(A*x_known-b)
                lastIter = j-1

                allVectors   = [obj.Yk V(:,(1:j))];
                totalVectors = size(allVectors,2);


                obj.Yk = allVectors(:,randperm(totalVectors,omega)); %%%%%!!!!!

                obj.x0 = x_k;
          
      
            end
        end
        


        function [c,s,S,Gold_1,Gnew] = expandS(obj,Tcol,Beta,Gold_1,Gnew,j)            
            Gold_2 = Gold_1; %G1
            Gold_1 = Gnew; %G2

            %expand G's

            Gold_2 = [Gold_2 zeros(j,1); zeros(1,j) 1]; %% fix this
            Gold_1 = [Gold_1 zeros(j,1); zeros(1,j) 1]; %% fix this

            Scol   = Gold_1*(Gold_2*Tcol);
            [c, s] = obj.computeGivensRotation(Scol);

            Gnew   = [eye(j-1) zeros(j-1,2);zeros(1,j-1) c s; zeros(1,j-1) s -c];

            S(1:j-1,j)  = Scol(1:j-1);
            S(j,j)      = c*Scol(j)+s*Beta;            
        end


        function [incX,incB,S,vTilde,bTilde,V,yHat,Alpha,Beta,Gold_1,Gnew] = Lanczos(obj,vHat,bCols,S,vTilde,bTilde,V,yHat,Alpha,Beta,Gold_1,Gnew,j)

            vHat        = vHat - V(:,j-1)*Beta(j);
            Alpha(j)    = V(:,j)'*vHat;
            vHat        = vHat - V(:,j)*Alpha(j);
            Beta(j+1)   = norm(vHat);
            V(:,j+1)    = vHat/Beta(j+1);


            Tcol = [zeros(j-2,1); Beta(j); Alpha(j); Beta(j+1)];
                        
            [c,s,S,Gold_1,Gnew] = obj.expandS(Tcol,Beta(j+1),Gold_1,Gnew,j);
            
            %S(j+1,j) = 0;

            modif       = yHat(j);
            yHat(j)     = c*modif;
            yHat(j+1)   = s*modif;

            vTilde(:,j) = (1/S(j,j))*(V(:,j)-vTilde(:,j-1)*S(j-1,j)-vTilde(:,j-2)*S(j-2,j)); %equivalent d
            bTilde(:,j) = (1/S(j,j))*(bCols(:,j)-bTilde(:,j-1)*S(j-1,j)-bTilde(:,j-2)*S(j-2,j));

            incX       = vTilde(:,j)*yHat(j);
            incB       = bTilde(:,j)*yHat(j);

        end

        function [c, s] = computeGivensRotation(obj, Scol)
            s1 = Scol(end-1);
            s2 = Scol(end);

            shi = s1^2 + s2^2 + sign(s1)*s1*sqrt(s1^2+s2^2);
            c   = ((s2^2)/shi)-1;
            s   = -(s2*(s1 + sign(s1)*sqrt(s1^2+s2^2)))/shi;
        end
    end

    methods (Access = private)


        function init(obj)       
        end

        function computeInitialValues(obj, A, b)
            if isempty(obj.x0)
                obj.x0 = zeros(obj.n,1);
            end
            obj.r0 = b-A*obj.x0;
        end

        function [Alpha, Beta]  = defineTVectors(obj)
            Alpha   = zeros(obj.n,1);
            Beta    = zeros(obj.n,1);
        end

        function [V, vTilde]    = defineVVectors(obj, k)
            V       = zeros(obj.n,k);
            vTilde  = zeros(obj.n,k);
        end

        function [bCols, bTilde]   = defineBVectors(obj, k, w)
            bTilde  = zeros(w,k);
            bCols   = zeros(w,k);
        end

        function [bHat, yHat] = defineRecyclingVectors(obj, Ck, R, res)
            bHat    = R\(Ck'*obj.r0);
            yHat    = zeros(obj.n,1);
            yHat(1) = res;
        end

        function Pk = createPk(obj, totalVectors, s)
            Pk = zeros(totalVectors,s);
            Pk(1:s,1:s) = eye(s);            
            Pk = Pk(randperm(totalVectors),:);
        end

        function  [c, s, Gamma_j] = computeH(Betajm1, Gamma_j)
            shi     = Gamma_j^2 + Betajm1^2 + sign(Gamma_j)*Gamma_j*sqrt(Gamma_j^2+Betajm1^2);
            c       = ((Betajm1^2)/shi)-1;
            s       = -(Betajm1*(Gamma_j + sign(Gamma_j)*sqrt(Gamma_j^2+Betajm1^2)))/shi;

            Gamma_j = c*Gamma_j + s*Betajm1;
        end



        function [Ck, R] = computeQR(obj, A);
            qrA      = A*obj.Yk;
            [Ck, R]  = qr(qrA, 0);
        end

    end
end

