classdef test2DSeqLaminateInVoigtWithFormuleOfAllaireWebPAgeExercise < test
    
    properties (Access = private)
        E
        nu
        lambda2D
        
        mu
        lambda
        
        ChAllaireWeb
        ChRank2
    end
    
    methods (Access = public)
        
        function obj = test2DSeqLaminateInVoigtWithFormuleOfAllaireWebPAgeExercise()
            obj.lambda = 0.7500;
            obj.mu = 0.3750;    
            obj.computeYoungFromMuAndLambda();
            obj.computePoissonFromMuAndLambda();
            obj.computeLambda2DFromYoungAndPoisson();
            
            d1 = [1 1 0];
            d2 = [rand(1,2)  0];
            d1 = d1/norm(d1);
            d2 = d2/norm(d2);
            m1 = 1;
            m2 = 0;
            theta = 0.3;
            epsil = 1e-3;
            
            obj.ChAllaireWeb = RankTwoLaminateHomogenizerFromAllaireTestedWebPageCode.ChinvOld(obj.lambda,obj.lambda2D,obj.mu,d1,d2,m1,m2,theta,epsil);
            Homog =  Rank2Homogenizer2D(obj.E,obj.nu,d1,d2,m1,m2,theta,epsil);
            obj.ChRank2 = Homog.Ch;
            
            Tens1 = IsotropicConstitutiveTensor3D(obj.E,obj.nu);            
            Tens1.computeTensorVoigtInPlaneStress();
            Tens0  = IsotropicConstitutiveTensor3D(obj.E*epsil,obj.nu);            
            Tens0.computeTensorVoigtInPlaneStress();
            
            R2 = RankTwoLaminateHomogenizer(Tens1,Tens0,[d1;d2],[m1 m2],theta);
            ChR2_2 = double(R2.Ch);
        end
        
        function computeYoungFromMuAndLambda(obj)
            obj.E = obj.mu*(3*obj.lambda+2*obj.mu)/(obj.lambda+obj.mu);
        end
        
        function computePoissonFromMuAndLambda(obj)
            obj.nu =  obj.lambda/(2*(obj.lambda + obj.mu));
        end
        
        function computeLambda2DFromYoungAndPoisson(obj)
            obj.lambda2D = obj.E*obj.nu/(1+obj.nu)/(1-obj.nu);
        end
        
    end
    
    methods (Access = protected)
            
        function hasPassed = hasPassed(obj)
            ChWeb = double(obj.ChAllaireWeb);
            ChR2 = double(obj.ChRank2);
            hasPassed = norm(ChWeb - ChR2)/norm(ChR2) < 1e-4;
        end
    end
end

