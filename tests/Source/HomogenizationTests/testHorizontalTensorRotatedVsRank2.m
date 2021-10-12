classdef testHorizontalTensorRotatedVsRank2 < ...
        testHorizontalTensorRotatedVsSequentialLaminate

    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsRank2()
            obj.computeTest()
        end
        
        function hasPassed = hasPassed(obj)
            rotHor = obj.rotHorTensor.getValue();
            lTens  = obj.lamTensor.getValue();
            hasPassed = norm(rotHor - lTens) > 1e-3;
        end

    end
    
    methods (Access = protected)

        function computeLaminateDirectly(obj)
            lambda2D = obj.C1.getLambda2D();
            mu    = obj.C1.getMu();
            c0    = obj.makeItPlaneStress(obj.C0);
            c1    = obj.makeItPlaneStress(obj.C1);
            dir{1}   = obj.lamDir;
            dir{2}   = obj.lamDir;
            mi(1)    = obj.lamPar;
            mi(2)    = 0;
            frac  = obj.theta;
            rank2 = RankTwoLaminateHomogenizer(c0,c1,dir,mi,frac,lambda2D,mu);
            obj.lamTensor = SymmetricFourthOrderPlaneStressVoigtTensor;
            obj.lamTensor.setValue(rank2.getTensor());
        end

    end
    
    methods (Access = private, Static)

        function c = makeItPlaneStress(c)
            c    = Tensor2VoigtConverter.convert(c);
            c    = PlaneStressTransformer.transform(c);
            c    = c.getValue();
        end

    end

end