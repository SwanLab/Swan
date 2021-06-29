classdef testHorizontalTensorRotatedVsVPH < ...
        testHorizontalTensorRotatedVsSequentialLaminate & ...
        testShowingError
    
    properties (Access = protected)
        tol = 1e-12;
    end
    
    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsVPH()
            obj.computeTest()
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
        
        function computeError(obj)
            rotHor = obj.rotHorTensor.getValue();
            lTens  = obj.lamTensor.getValue();
            obj.error = norm(rotHor - lTens);
        end
        
    end
    
end

