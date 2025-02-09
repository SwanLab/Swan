classdef RankTwoLaminateHomogenizerFromAllaireTestedWebPageCode
    
    properties
    end
    
    methods (Static)
        function Chomog = ChinvOld(lambda,mu,d1,d2,m1,m2,theta,epsil)
            
            sqrt2 = sqrt(2);
            
            e1x = d1(1);
            e1y = d1(2);
            
            e2x = d2(1);
            e2y = d2(2);

            
            K = (mu+lambda)/(mu*(2*mu+lambda));
            
            muVoigt = mu;
            
            A =	m1 * ( (lambda+2*mu) - 1/mu*(lambda^2*e1y^2+(lambda+2*mu)^2*e1x^2) + K*((lambda+2*mu)*e1x^2+lambda*e1y^2)^2 ) + ...
                +m2* ( (lambda+2*mu) - 1/mu*(lambda^2*e2y^2+(lambda+2*mu)^2*e2x^2) + K*((lambda+2*mu)*e2x^2+lambda*e2y^2)^2 );
            B =	m1 * ( (lambda+2*mu) - 1/mu*(lambda^2*e1x^2+(lambda+2*mu)^2*e1y^2) + K*((lambda+2*mu)*e1y^2+lambda*e1x^2)^2 )  + ...
                +m2* ( (lambda+2*mu) - 1/mu*(lambda^2*e2x^2+(lambda+2*mu)^2*e2y^2) + K*((lambda+2*mu)*e2y^2+lambda*e2x^2)^2 );
            C =	m1 * ( 4*mu - 1/mu*(2*mu)^2 + K*(4*mu*e1x*e1y)^2 )/2  + ...
                +m2* ( 4*mu - 1/mu*(2*mu)^2 + K*(4*mu*e2x*e2y)^2 )/2;
            D = 	m1 * ( 2*lambda - 1/mu*(2*lambda*(lambda+2*mu)) + K*2*((lambda+2*mu)*e1y^2+lambda*e1x^2)*((lambda+2*mu)*e1x^2+lambda*e1y^2) )/2  + ...
                +m2* ( 2*lambda - 1/mu*(2*lambda*(lambda+2*mu)) + K*2*((lambda+2*mu)*e2y^2+lambda*e2x^2)*((lambda+2*mu)*e2x^2+lambda*e2y^2) )/2;
            E =	m1 * ( -1/mu*(4*mu*(2*lambda+2*mu)*e1x*e1y) + K*2*(4*mu*e1x*e1y*((lambda+2*mu)*e1y^2+lambda*e1x^2)) )/(2*sqrt2)  + ...
                +m2* ( -1/mu*(4*mu*(2*lambda+2*mu)*e2x*e2y) + K*2*(4*mu*e2x*e2y*((lambda+2*mu)*e2y^2+lambda*e2x^2)) )/(2*sqrt2);
            F = 	m1 * ( -1/mu*(4*mu*(2*lambda+2*mu)*e1x*e1y) + K*2*(4*mu*e1x*e1y*((lambda+2*mu)*e1x^2+lambda*e1y^2)) )/(2*sqrt2)  + ...
                +m2* ( -1/mu*(4*mu*(2*lambda+2*mu)*e2x*e2y) + K*2*(4*mu*e2x*e2y*((lambda+2*mu)*e2x^2+lambda*e2y^2)) )/(2*sqrt2);
            
            %// Ajout du matériau mou (qui simule le vide)
            A = epsil/(1.-epsil)*(lambda+2*mu)	+ theta*A;
            B = epsil/(1.-epsil)*(lambda+2*mu)	+ theta*B;
            C = epsil/(1.-epsil)*muVoigt		+ theta*C;
            D = epsil/(1.-epsil)*lambda		+ theta*D;
            E = epsil/(1.-epsil)*0			+ theta*E;
            F = epsil/(1.-epsil)*0			+ theta*F;
            
            %//Première inversion
            DET = A*B*C-A*E*E-B*F*F-C*D*D+2*D*E*F;
            A1=(B*C-E*E)/DET;
            B1=(A*C-F*F)/DET;
            C1=(A*B-D*D)/DET;
            D1=(E*F-C*D)/DET;
            E1=(D*F-A*E)/DET;
            F1=(D*E-B*F)/DET;
            
            A = (1-theta)*A1+(lambda+2*mu)/(4*mu*(lambda+mu));
            B = (1-theta)*B1+(lambda+2*mu)/(4*mu*(lambda+mu));
            C = (1-theta)*C1+1/(muVoigt);
            D = (1-theta)*D1-lambda/(4*mu*(lambda+mu));
            E = (1-theta)*E1;
            F = (1-theta)*F1;
            
            %//Deuxième inversion
            DET = A*B*C-A*E*E-B*F*F-C*D*D+2*D*E*F;
            A1=(B*C-E*E)/DET;
            B1=(A*C-F*F)/DET;
            C1=(A*B-D*D)/DET;
            D1=(E*F-C*D)/DET;
            E1=(D*F-A*E)/DET;
            F1=(D*E-B*F)/DET;
            A=A1;
            B=B1;
            C=C1;
            D=D1;
            E=E1;
            F=F1;
            
            Chomog = [A   D	 F;
                D   B  E;
                F   E  C];
            
        end
        
    end
    
end