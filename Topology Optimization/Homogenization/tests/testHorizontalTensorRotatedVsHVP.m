classdef testHorizontalTensorRotatedVsHVP < ...
        testHorizontalTensorRotatedVsSequentialLaminate
    
    properties (Access = protected)
        tol
    end
    
    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsHVP()
            obj.computeTest()
        end
    end
    
    methods (Access = protected)
        
        function computeLaminateDirectly(obj)
            c0       = obj.C0;
            c1       = obj.C1;
            dir(1,:) = obj.lamDir;
            m1       = obj.lamPar;
            frac     = obj.theta;
            lam      = HomogVoigtPlaneStressHomogenizer(c0,c1,dir,m1,frac);
            obj.lamTensor = lam.getPlaneStressHomogenizedTensor();
        end
        
        function createTolerance(obj)
           obj.tol = 1e-3; 
        end
    end
    
end

