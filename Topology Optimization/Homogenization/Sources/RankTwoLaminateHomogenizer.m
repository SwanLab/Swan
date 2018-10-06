classdef RankTwoLaminateHomogenizer
    
    properties (Access = private)
        FirstDirection
        SecondDirection
        Fraction
        StiffTensor
        WeakTensor
        Kparameter
        AnisotropicTensor
        
        mu
        lambda2D
    end
    
    properties (Access = public)
        Ch
    end
    
    methods 
        
       function obj = RankTwoLaminateHomogenizer(StiffTensor,WeakTensor,Directions,LamParams,theta)
            
            d1 = Directions(1,:);
            d2 = Directions(2,:);
            
            m1 = LamParams(1);
            m2 = LamParams(2);
            
            obj.Fraction = theta;
            obj.StiffTensor = StiffTensor;
            obj.WeakTensor = WeakTensor;
            
            
            
            C1 = obj.StiffTensor.tensorVoigtInPlaneStress; 

            C0 = obj.WeakTensor.tensorVoigtInPlaneStress;
           
           
            S1 = obj.InvertSymmMatrix(C1);
            
            
            S0 = obj.InvertSymmMatrix(C0);
            
            obj.mu = StiffTensor.mu;
            
            
            lambda2D = obj.StiffTensor.E*obj.StiffTensor.nu/(1+obj.StiffTensor.nu)/(1-obj.StiffTensor.nu);
            obj.lambda2D = lambda2D;
            %obj.lambda =  StiffTensor.lambda;    
            
            
            obj.Kparameter = (obj.mu+obj.lambda2D)/(obj.mu*(2*obj.mu+obj.lambda2D));
            Cm1 = obj.computeChCorrector(d1);
            Cm2 = obj.computeChCorrector(d2);
            
            obj.AnisotropicTensor{1} = fourthOrderTensor();
            obj.AnisotropicTensor{2} = fourthOrderTensor();
            
            
            obj.AnisotropicTensor{1}.tensorVoigtInPlaneStress = Cm1;
            obj.AnisotropicTensor{2}.tensorVoigtInPlaneStress = Cm2;
            
            Cm = Cm1*m1 + Cm2*m2;
            
            S01 = S0 - S1;
            C01 = obj.InvertSymmMatrix(S01);
            
            Ctheta = C01 + theta*Cm;   
            
            Stheta = obj.InvertSymmMatrix(Ctheta);
            
            Sh = S1 +(1-theta)*Stheta;            
            obj.Ch = obj.InvertSymmMatrix(Sh);
            
         end 
        
        
        function val = computeOneOfTheTwoFirstDiagonalTerms(obj,ex,ey)
            lam = obj.lambda2D;
            muT = obj.mu;
            K = obj.Kparameter;
            val = ( (lam+2*muT) - 1/muT*(lam^2*ey^2+(lam+2*muT)^2*ex^2) + K*((lam+2*muT)*ex^2+lam*ey^2)^2 );            
        end
        
        function val = computeCrossThirdTerms(obj,ex,ey)
            lam = obj.lambda2D;
            muT = obj.StiffTensor.mu;
            K = obj.Kparameter;
            val = ( -1/muT*(4*muT*(2*lam+2*muT)*ex*ey) + K*2*(4*muT*ex*ey*((lam+2*muT)*ey^2+lam*ex^2)) )/(2*sqrt(2));
        end
        
        function C11 = computeC11(obj,direction)
            ex = direction(1);
            ey = direction(2);
            C11 = obj.computeOneOfTheTwoFirstDiagonalTerms(ex,ey);
        end

        function C22 = computeC22(obj,direction)
            ex = direction(1);
            ey = direction(2);
            C22 = obj.computeOneOfTheTwoFirstDiagonalTerms(ey,ex);
        end
        
        function C33 = computeC33(obj,direction)
            ex = direction(1);
            ey = direction(2);
            muT = obj.mu;
            K = obj.Kparameter;
            C33 = ( 4*muT - 1/muT*(2*muT)^2 + K*(4*muT*ex*ey)^2 )/2;            
        end
        

        function C21 = computeC21(obj,direction) 
            ex = direction(1);
            ey = direction(2);
            lam = obj.lambda2D;
            muT = obj.mu;
            K = obj.Kparameter;
            C21 = ( 2*lam - 1/muT*(2*lam*(lam+2*muT)) + K*2*((lam+2*muT)*ey^2+lam*ex^2)*((lam+2*muT)*ex^2+lam*ey^2) )/2;
        end
        
        function C23 = computeC23(obj,direction) 
            ex = direction(1);
            ey = direction(2);
            C23 = obj.computeCrossThirdTerms(ex,ey);
        end
        
        function C13 = computeC13(obj,direction) 
            ex = direction(1);
            ey = direction(2);
            C13 = obj.computeCrossThirdTerms(ey,ex);
        end
        
        function Chomog = computeChCorrector(obj,direction)
            C11 = obj.computeC11(direction);
            C22 = obj.computeC22(direction);
            C33 = obj.computeC33(direction);
            C21 = obj.computeC21(direction);
            C23 = obj.computeC23(direction);
            C13 = obj.computeC13(direction);
            
             Chomog = [ C11  C21 C13; 
                        C21  C22 C23;
                        C13  C23 C33];
        end
    end
    
    methods (Static)
         
        function InvCh = InvertSymmMatrix(Ch)
            C11 = Ch(1,1);
            C22 = Ch(2,2);
            C33 = Ch(3,3);
            C21 = Ch(2,1);
            C23 = Ch(2,3);
            C13 = Ch(1,3);
            
            DET = C11*C22*C33-C11*C23*C23-C22*C13*C13-C33*C21*C21+2*C21*C23*C13;
            InvC11 = (C22*C33-C23*C23)/DET;
            InvC22 = (C11*C33-C13*C13)/DET;
            InvC33 = (C11*C22-C21*C21)/DET;
            InvC21 = (C23*C13-C33*C21)/DET;
            InvC23 = (C21*C13-C11*C23)/DET;
            InvC13 = (C21*C23-C22*C13)/DET;
            
            InvCh = [InvC11   InvC21  InvC13;
                     InvC21   InvC22  InvC23;
                     InvC13   InvC23  InvC33];
            
        end
         
          
        
    end
    
end

