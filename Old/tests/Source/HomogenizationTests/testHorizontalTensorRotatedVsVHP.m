classdef testHorizontalTensorRotatedVsVHP < ...
        testHorizontalTensorRotatedVsSequentialLaminate

    methods (Access = public)

        function obj = testHorizontalTensorRotatedVsVHP()
            obj.computeTest()
        end

        function hasPassed = hasPassed(obj)
            rotHor = obj.rotHorTensor.getValue();
            lTens  = obj.lamTensor.getValue();
            lb = norm(rotHor - lTens) > 1e-12;
            ub = norm(rotHor - lTens) < 1e-2;
            hasPassed = lb & ub;
        end

    end
    
    methods (Access = protected)

        function computeLaminateDirectly(obj)
            c0       = obj.C0;
            c1       = obj.C1;
            dir{1}   = obj.lamDir;
            m1       = obj.lamPar;
            frac     = obj.theta;
            lam      = VoigtHomogPlaneStressHomogenizer(c0,c1,dir,m1,frac);
            obj.lamTensor = lam.getPlaneStressHomogenizedTensor();
        end

    end

end