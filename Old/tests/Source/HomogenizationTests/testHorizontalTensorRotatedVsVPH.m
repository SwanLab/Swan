classdef testHorizontalTensorRotatedVsVPH < ...
        testHorizontalTensorRotatedVsSequentialLaminate
    
    properties (Access = public)
        tol = 1e-12;
    end
    
    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsVPH()
            obj.computeTest();
        end

        function error = computeError(obj)
            rotHor = obj.rotHorTensor.getValue();
            lTens  = obj.lamTensor.getValue();
            error = norm(rotHor - lTens);
        end
    end
    
    methods (Access = protected)
        
        function computeLaminateDirectly(obj)
            c0       = obj.C0;
            c1       = obj.C1;
            dir{1}   = obj.lamDir;
            m1       = obj.lamPar;
            frac     = obj.theta;
            lam      = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,m1,frac);
            obj.lamTensor = lam.getPlaneStressHomogenizedTensor();
        end

    end

end