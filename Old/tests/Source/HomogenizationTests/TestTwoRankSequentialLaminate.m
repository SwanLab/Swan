classdef TestTwoRankSequentialLaminate < handle

    properties (Access = public)
        tol = 1e-12;
    end

    properties (Access = private)
      C0
      C1
      mu
      lambda
      lambda2D
      FractionVolume
      d1
      d2
      m1
      m2
      FirstCheckedRank2Ch
      WebRank2Ch
      Rank2Ch
    end
    
    methods (Access = public)
        
        function obj = TestTwoRankSequentialLaminate()
            obj.createTensors();
            obj.createParameters();
            obj.loadHomogenizedCheckedConstitutiveTensor();
            obj.computeTwoRankSequentialLaminate();
            obj.computeTwoRankSequentialLaminateFromAllaireTestedWebPageCode();
        end

        function error = computeError(obj)
            InitCh  = double(obj.FirstCheckedRank2Ch);
            R2Ch    = double(obj.Rank2Ch);
            R2ChWeb = double(obj.WebRank2Ch);
            err1  = norm(InitCh - R2Ch)/norm(R2Ch);
            err2  = norm(R2ChWeb- R2Ch)/norm(R2Ch);
            error = max(err1,err2);
        end

    end

    methods (Access = private)

        function createTensors(obj)
            obj.createStiffTensor();
            obj.createWeakTensor();
            obj.storeLambda2D();
            obj.makeTensorsVoigtPlaneStress();
        end
        
        function createStiffTensor(obj)
            obj.lambda = 0.7500;
            obj.mu     = 0.3750;
            obj.C1 = IsotropicConstitutiveTensor.createWithLambdaAndMu(obj.lambda,obj.mu);
        end
        
        function createWeakTensor(obj)
            E  = obj.C1.getYoung();
            nu = obj.C1.getPoisson();
            epsil = 1.0000e-03;
            E0 = epsil*E;
            nu0 = nu;
            obj.C0 = IsotropicConstitutiveTensor(E0,nu0);
        end
        
        function storeLambda2D(obj)
            obj.lambda2D = obj.C1.getLambda2D();
        end

        function makeTensorsVoigtPlaneStress(obj)
            obj.C0 = obj.makeTensorVoigtPlaneStress(obj.C0);
            obj.C1 = obj.makeTensorVoigtPlaneStress(obj.C1);
        end
        
        function createParameters(obj)
            obj.loadDirections()
            obj.loadLaminationParameters()
            obj.loadFractionVolume()
        end
        
        function loadDirections(obj)
            obj.loadFirstDirection()
            obj.loadSecondDirection()
        end
        
        function loadFirstDirection(obj)
            dir = [1     0     0];
            obj.d1 = Vector3D;
            obj.d1.setValue(dir);
            obj.d1.normalize()
        end
        
        function loadSecondDirection(obj)
            dir = [1     3     2];
            obj.d2 = Vector3D;
            obj.d2.setValue(dir);
            obj.d2.normalize()
        end
        
        function loadLaminationParameters(obj)
            obj.m1  = 0.8;
            obj.m2  = 0.2;
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
            dir{1} = obj.d1;
            dir{2} = obj.d2;
            params = [obj.m1;obj.m2];
            theta = obj.FractionVolume;
            c1    = obj.C1.getValue();
            c0    = obj.C0.getValue();
            lam2D = obj.lambda2D;
            muV   = obj.mu;
            Homogenizer = RankTwoLaminateHomogenizer(c1,c0,dir,params,...
                                                     theta,lam2D,muV);
            obj.Rank2Ch = Homogenizer.getTensor();
        end
        
        
        function computeTwoRankSequentialLaminateFromAllaireTestedWebPageCode(obj)
            %Lambda   = obj.lambda2D;
            Lambda2D = obj.lambda2D;
            muV = obj.mu();
            theta = obj.FractionVolume;
            epsil = 1.0000e-03;
            dir1 = obj.d1.getValue();
            dir2 = obj.d2.getValue();
            obj.WebRank2Ch = RankTwoLaminateHomogenizerFromAllaireTestedWebPageCode.ChinvOld(...
                             Lambda2D,muV,dir1,dir2,obj.m1,obj.m2,theta,epsil);
        end
        
    end
    
    methods (Access = private, Static)
               
        function CVoigtPS = makeTensorVoigtPlaneStress(C)
            CVoigt = Tensor2VoigtConverter.convert(C);
            CVoigtPS = PlaneStressTransformer.transform(CVoigt);
        end
        
    end

end