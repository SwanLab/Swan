classdef TestGeneralTwoRankSequentialLaminate < handle
    
    properties (Access = private)
      fractionVolume
      directions
      lamParams
      stiffTensor
      weakTensor
      Rank2Ch
      SeqLamCh
    end
    
    methods (Access = public)

        function obj = TestGeneralTwoRankSequentialLaminate()
            obj.init()
            obj.computeTwoRankSequentialLaminate()
            obj.computeGeneralTwoRankSequentialLaminate()
        end

        function hasPassed = hasPassed(obj)
            RankTwoCh  = obj.Rank2Ch;
            SqCh = obj.SeqLamCh.getValue();
            hasPassed = norm(RankTwoCh - SqCh)/norm(SqCh) > 1e-6;
        end

    end

    methods (Access = private)

        function init(obj)
            obj.fractionVolume = 0.8000;
            obj.lamParams = [1 0];
            obj.loadLaminateDirections()
            obj.createStiffAndWeakTensors()
        end

        function loadLaminateDirections(obj)
            d1 = [rand(1) rand(1)  0];
            d2 = [1     3     0];
            obj.directions{1} = obj.createDirection(d1);
            obj.directions{2} = obj.createDirection(d2);
        end

        function createStiffAndWeakTensors(obj)
            epsilon = 1.0000e-03;
            E1  = 1;
            nu1 = 1/3;
            E0  = epsilon*E1;
            nu0 = 1/3;
            obj.stiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.weakTensor  = IsotropicConstitutiveTensor(E0,nu0);
        end

        function computeTwoRankSequentialLaminate(obj)
            mi    = obj.lamParams;
            dir   = obj.directions;
            theta = obj.fractionVolume;
            C1    = obj.stiffTensor.getValue();
            C0    = obj.weakTensor.getValue();
            mu    = obj.stiffTensor.getMu();
            lambda2D = obj.stiffTensor.getLambda2D();
            homogenizer = RankTwoLaminateHomogenizer(C1,C0,dir,mi,theta,lambda2D,mu);
            obj.Rank2Ch = homogenizer.getTensor;
        end
        
        function computeGeneralTwoRankSequentialLaminate(obj)
            C0       = obj.weakTensor;
            C1       = obj.stiffTensor;
            mi       = obj.lamParams;
            dir      = obj.directions;
            Theta    = obj.fractionVolume;
            SeqHomog      = VoigtPlaneStressHomogHomogenizer(C0,C1,dir,mi,Theta);
            obj.SeqLamCh  = SeqHomog.getPlaneStressHomogenizedTensor();
        end

    end

    methods (Access = private, Static)

        function dir = createDirection(d)
            dir = Vector3D;
            dir.setValue(d);
            dir.normalize();
        end

    end
 
end