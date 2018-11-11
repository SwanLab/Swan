classdef TestTwoRankSequentialLaminate < test
    
    
    
    properties (Access = private)
      lambda
      mu
      E 
      nu
      lambda2D
      
      FractionVolume      
      
      FirstLaminateDirection
      SecondLaminateDirection
      
      FirstLaminationParameter
      SecondLaminationParameter
      
      FirstCheckedRank2Ch
      WebRank2Ch
      Rank2Ch

    end
    
    methods (Access = public) 
        
        function obj = TestTwoRankSequentialLaminate()
            obj.loadLameParameters()
            obj.loadDirections()
            obj.loadLaminationParameters()
            obj.loadFractionVolume()
            obj.loadHomogenizedCheckedConstitutiveTensor()
            obj.computeTwoRankSequentialLaminate()
            obj.computeTwoRankSequentialLaminateFromAllaireTestedWebPageCode()
        end
        
        
    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            InitCh  = double(obj.FirstCheckedRank2Ch);
            R2Ch    = double(obj.Rank2Ch);
            R2ChWeb = double(obj.WebRank2Ch); 
            
            firstCheck  = norm(InitCh - R2Ch )/norm(R2Ch) < 1e-4;
            secondCheck = norm(R2ChWeb- R2Ch)/norm(R2Ch) < 1e-4;
            hasPassed = firstCheck && secondCheck;
        end
    end
    
    methods (Access = private)
        
        function loadLameParameters(obj)
            obj.lambda = 0.7500;
            obj.mu = 0.3750;     
            obj.computeYoungFromMuAndLambda();
            obj.computePoissonFromMuAndLambda();
            obj.computeLambda2DFromYoungAndPoisson();
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
        
        function loadDirections(obj)
            obj.loadFirstDirection()
            obj.loadSecondDirection()
        end
        
        function loadFirstDirection(obj)
            direction = [1     0     0];
            obj.FirstLaminateDirection = obj.normalizeDirection(direction);
        end
        
        function loadSecondDirection(obj)
            direction = [1     3     2];
            obj.SecondLaminateDirection = obj.normalizeDirection(direction);
        end
        
        function loadLaminationParameters(obj)
            obj.FirstLaminationParameter  = 0.8;
            obj.SecondLaminationParameter = 0.2;
        end
        
        function loadFractionVolume(obj)
            obj.FractionVolume = 0.8000;
        end
        
        function loadHomogenizedCheckedConstitutiveTensor(obj)
            obj.FirstCheckedRank2Ch = [    0.326824930177257   0.028535636891760  -0.155777642710044
                                           0.028535636891760   0.746824571833621  -0.074187675387998
                                          -0.155777642710044  -0.074187675387998   0.032148130248805];
        end
        
        function computeTwoRankSequentialLaminate(obj)           
            d1 = obj.FirstLaminateDirection;
            d2 = obj.SecondLaminateDirection;
            m1 = obj.FirstLaminationParameter;
            m2 = obj.SecondLaminationParameter;
            
            Dir    = [d1;d2];
            Params = [m1 m2];
            
            Theta = obj.FractionVolume;
            
            epsil = 1.0000e-03;
            E1 = obj.E;
            nu1 = obj.nu;
            E0 = epsil*obj.E;
            nu0 = obj.nu;
            
            C1 = IsotropicConstitutiveTensor3D(E1,nu1);            
            C1.computeTensorVoigtInPlaneStress()
            C0  = IsotropicConstitutiveTensor3D(E0,nu0);            
            C0.computeTensorVoigtInPlaneStress()

            
            Homogenizer = RankTwoLaminateHomogenizer(C1,C0,Dir,Params,Theta);
            obj.Rank2Ch = Homogenizer.Ch;
        end
        
        function computeTwoRankSequentialLaminateFromAllaireTestedWebPageCode(obj)
            Lambda = obj.lambda;
            Lambda2D = obj.lambda2D;
            Mu = obj.mu;
            d1 = obj.FirstLaminateDirection;
            d2 = obj.SecondLaminateDirection;
            m1 = obj.FirstLaminationParameter;
            m2 = obj.SecondLaminationParameter;
            theta = obj.FractionVolume;
            epsil = 1.0000e-03;
            obj.WebRank2Ch = RankTwoLaminateHomogenizerFromAllaireTestedWebPageCode.ChinvOld(Lambda2D,Mu,d1,d2,m1,m2,theta,epsil);
        end
        
    end
    
    methods (Access = private, Static)
        function direction = normalizeDirection(direction)
            direction = direction/norm(direction);
        end
    end
    
    
    
end