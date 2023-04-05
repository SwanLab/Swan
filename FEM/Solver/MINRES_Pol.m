classdef MINRES_Pol < Solver


        
      methods (Static)

        function x_k = solve(A,b)
        n = size(A,1);
        k = 3000;
        
        Alpha = zeros(k+1,1);
        
        Beta    = zeros(k+1,1);
        Beta(1) = norm(b);
        
        V       = zeros(n,k);
        v_0     = 0;
        V(:,1)  = b/Beta(1);
        
        c = -1;
        s = 0;
        
        Gamma   = zeros(k,1);
        Delta   = zeros(k,1);
        Epsilon = zeros(k,1);
        
        Gamma(1) = 0;
        
        phi = Beta(1);
        
        d = zeros(n,k);
        x0   = 0;
        
        t = zeros(k,1);
        
        
        %j=1
            j = 1;
            x_prev = x0;
            [Alpha(j), Beta(j+1), V(:,j+1)] = lanczosProcess3(A, V(:,j), v_0, Beta(j));
            
            deltaj_prev = Delta(j);
            Delta(j) = c*deltaj_prev + s*Alpha(j);
            Gamma(j) = s*deltaj_prev - c*Alpha(j);
            Epsilon(j+1) = s*Beta(j+1);
            Delta(j+1) = -c*Beta(j+1);
        
            shi = Gamma(j)^2 + Beta(j+1)^2 + sign(Gamma(j))*Gamma(j)*sqrt(Gamma(j)^2+Beta(j+1)^2);
            c = ((Beta(j+1)^2)/shi)-1;
            s = -(Beta(j+1)*(Gamma(j) + sign(Gamma(j))*sqrt(Gamma(j)^2+Beta(j+1)^2)))/shi;
            
            phi_prev = phi;
            phi  = s*phi_prev;
            t(j) = c*phi_prev;
        
            Gamma(j) = c*Gamma(j) + s*Beta(j+1);
            % R(:,j) is now completed
        
            d(:,j) = V(:,j)/Gamma(j);
        
            x1 = x_prev + t(j)*d(:,j);
        
        %j=2
            j = 2;
            
            x_prev = x1;
        
            [Alpha(j), Beta(j+1), V(:,j+1)] = lanczosProcess3(A, V(:,j), V(:,j-1), Beta(j));
            
            deltaj_prev = Delta(j);
            Delta(j) = c*deltaj_prev + s*Alpha(j);
            Gamma(j) = s*deltaj_prev - c*Alpha(j);
            Epsilon(j+1) = s*Beta(j+1);
            Delta(j+1) = -c*Beta(j+1);
        
            shi = Gamma(j)^2 + Beta(j+1)^2 + sign(Gamma(j))*Gamma(j)*sqrt(Gamma(j)^2+Beta(j+1)^2);
            c = ((Beta(j+1)^2)/shi)-1;
            s = -(Beta(j+1)*(Gamma(j) + sign(Gamma(j))*sqrt(Gamma(j)^2+Beta(j+1)^2)))/shi;
            
            phi_prev = phi;
            phi  = s*phi_prev;
            t(j) = c*phi_prev;
        
            Gamma(j) = c*Gamma(j) + s*Beta(j+1);
            % R(:,j) is now completed
        
            d(:,j) = (V(:,j)-Delta(j)*d(:,j-1))/Gamma(j);
        
            xj = x_prev + t(j)*d(:,j);
        
            for j = 3:k
        
            x_prev = xj;
        
            [Alpha(j), Beta(j+1), V(:,j+1)] = lanczosProcess3(A, V(:,j), V(:,j-1), Beta(j));
            
            deltaj_prev = Delta(j);
            Delta(j) = c*deltaj_prev + s*Alpha(j);
            Gamma(j) = s*deltaj_prev - c*Alpha(j);
            Epsilon(j+1) = s*Beta(j+1);
            Delta(j+1) = -c*Beta(j+1);
        
            shi = Gamma(j)^2 + Beta(j+1)^2 + sign(Gamma(j))*Gamma(j)*sqrt(Gamma(j)^2+Beta(j+1)^2);
            c = ((Beta(j+1)^2)/shi)-1;
            s = -(Beta(j+1)*(Gamma(j) + sign(Gamma(j))*sqrt(Gamma(j)^2+Beta(j+1)^2)))/shi;
            
            phi_prev = phi;
            phi  = s*phi_prev;
            t(j) = c*phi_prev;
        
            Gamma(j) = c*Gamma(j) + s*Beta(j+1);
            % R(:,j) is now completed
        
            d(:,j) = (V(:,j)-Delta(j)*d(:,j-1)-Epsilon(j)*d(:,j-2))/Gamma(j);
        
            xj = x_prev + t(j)*d(:,j); 
            end
        
            x_k = xj;

            res = norm(A*x_k-b)/norm(b);
        end
      end
end