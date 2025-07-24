classdef AnisotropicContributionTensorForRank2 < handle
    
    properties (Access = private)
        dir
        mu
        lambda2D
        Kparameter
        Cani
    end
    
    methods (Access = public)
        
        function obj = AnisotropicContributionTensorForRank2(dir,mu,lambda2D)
            obj.init(dir,mu,lambda2D)
            obj.computeKparameter()
            obj.computeAnisotropicContributionTensor()
        end
        
        function C = getTensor(obj)
            C = obj.Cani;
        end
    end
    
    methods (Access = private)
        
        function init(obj,dir,mu,lambda2D)
            obj.dir = dir;
            obj.mu = mu;
            obj.lambda2D = lambda2D;
        end
        
        function computeKparameter(obj)
            muV = obj.mu;
            lam2D = obj.lambda2D;
            obj.Kparameter =  (muV+lam2D)/(muV*(2*muV+lam2D));
        end
        
        function computeAnisotropicContributionTensor(obj)
            C11 = obj.computeC11();
            C22 = obj.computeC22();
            C33 = obj.computeC33();
            C21 = obj.computeC21();
            C23 = obj.computeC23();
            C13 = obj.computeC13();
            
            Chomog = [ C11  C21 C13;
                       C21  C22 C23;
                       C13  C23 C33];
            obj.Cani = Chomog;
        end
        
         function C11 = computeC11(obj)
            ex = obj.dir(1);
            ey = obj.dir(2);
            C11 = obj.computeOneOfTheTwoFirstDiagonalTerms(ex,ey);
        end

        function C22 = computeC22(obj)
            ex = obj.dir(1);
            ey = obj.dir(2);
            C22 = obj.computeOneOfTheTwoFirstDiagonalTerms(ey,ex);
        end
        
        function C33 = computeC33(obj)
            ex = obj.dir(1);
            ey = obj.dir(2);
            muT = obj.mu;
            K = obj.Kparameter;
            C33 = ( 4*muT - 1/muT*(2*muT)^2 + K*(4*muT*ex*ey)^2 )/2;
        end
        
        function C21 = computeC21(obj) 
            ex = obj.dir(1);
            ey = obj.dir(2);
            lam = obj.lambda2D;
            muT = obj.mu;
            K = obj.Kparameter;
            C21 = ( 2*lam - 1/muT*(2*lam*(lam+2*muT)) + K*2*((lam+2*muT)*ey^2+lam*ex^2)*((lam+2*muT)*ex^2+lam*ey^2) )/2;
        end
        
        function C23 = computeC23(obj)
            ex = obj.dir(1);
            ey = obj.dir(2);
            C23 = obj.computeCrossThirdTerms(ex,ey);
        end
        
        function C13 = computeC13(obj) 
            ex = obj.dir(1);
            ey = obj.dir(2);
            C13 = obj.computeCrossThirdTerms(ey,ex);
        end
        
        function val = computeOneOfTheTwoFirstDiagonalTerms(obj,ex,ey)
            lam = obj.lambda2D;
            muT = obj.mu;
            K = obj.Kparameter;
            val = ( (lam+2*muT) - 1/muT*(lam^2*ey^2+(lam+2*muT)^2*ex^2) + K*((lam+2*muT)*ex^2+lam*ey^2)^2 );
        end
        
        function val = computeCrossThirdTerms(obj,ex,ey)
            lam = obj.lambda2D;
            muT = obj.mu;
            K = obj.Kparameter;
            val = ( -1/muT*(4*muT*(2*lam+2*muT)*ex*ey) + K*2*(4*muT*ex*ey*((lam+2*muT)*ey^2+lam*ex^2)) )/(2*sqrt(2));
        end
        
    end
    
end

