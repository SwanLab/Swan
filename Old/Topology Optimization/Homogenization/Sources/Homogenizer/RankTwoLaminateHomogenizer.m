classdef RankTwoLaminateHomogenizer < handle
    
    properties (Access = private)
        ani1
        ani2
        mi
        theta
        dir

        mu
        lambda2D
        
        C1
        C0
        homogTensor
        aniTensor
        incremTensor
        
        invMatTensor
        invFibTensor
        invThetaTensor
    end
    

    
    methods (Access = public)
        
       function obj = RankTwoLaminateHomogenizer(C1,C0,dir,mi,theta,lambda2D,mu)
            obj.init(C1,C0,dir,mi,theta,lambda2D,mu)
            obj.computeMatrixAndFiberComplianceTensors()
            obj.computeAnisotropicTensors()
            obj.computeIncrementalTensor()
            obj.computeAnisotropicContributionTensor();
            obj.computeThetaTensor()
            obj.computeHomogenizedTensor()
       end 
       
       function C = getTensor(obj)
           C = obj.homogTensor;
       end
       
    end
    
    methods (Access = private)
        
       function init(obj,C1,C0,dir,mi,theta,lambda2D,mu)
           obj.C1 = C1;
           obj.C0 = C0;
           obj.dir = dir;
           obj.mi = mi;
           obj.theta = theta;
           obj.lambda2D = lambda2D;
           obj.mu = mu;
       end
       
        function computeMatrixAndFiberComplianceTensors(obj)
            obj.invMatTensor = obj.invertTensor(obj.C0);
            obj.invFibTensor = obj.invertTensor(obj.C1);
        end
       
       function computeIncrementalTensor(obj)
            S0  = obj.invMatTensor;
            S1  = obj.invFibTensor;
            S01 = S0 - S1;
            C01 = obj.invertTensor(S01);
            obj.incremTensor = C01;
       end
        
      function computeAnisotropicTensors(obj)
            d1 = obj.dir{1};
            d2 = obj.dir{2};
            obj.ani1 = obj.computeAnisotropicTensor(d1);
            obj.ani2 = obj.computeAnisotropicTensor(d2);
      end
         
      function computeAnisotropicContributionTensor(obj)
            m1 = obj.mi(1);
            m2 = obj.mi(2);
            Cm1 = obj.ani1;
            Cm2 = obj.ani2;
            Cm = Cm1*m1 + Cm2*m2;
            obj.aniTensor = Cm;
      end
       
      function Cm = computeAnisotropicTensor(obj,dir)
          muV   = obj.mu;
          lam2D = obj.lambda2D;
          d = dir.getValue();
          aCont = AnisotropicContributionTensorForRank2(d,muV,lam2D);
          Cm = aCont.getTensor();
      end
      
      function computeThetaTensor(obj)
         C01 = obj.incremTensor;
         Ca  = obj.aniTensor;
         Ctheta = C01 + obj.theta*Ca;
         Stheta = obj.invertTensor(Ctheta);
         obj.invThetaTensor = Stheta;
      end
      
      function computeHomogenizedTensor(obj)
          S1     = obj.invFibTensor;
          Stheta = obj.invThetaTensor;
          Sh = S1 + (1-obj.theta)*Stheta;
          Ch = obj.invertTensor(Sh);
          obj.homogTensor = Ch;
      end

    end
    
    methods (Static)
         
        function InvCh = invertTensor(Ch)
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

