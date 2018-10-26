classdef testHorizontalTensorRotatedVsVHP < ...
        testHorizontalTensorRotatedVsSequentialLaminate
    
    properties (Access = protected)
        tol
    end
    
    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsVHP()
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
            lam      = VoigtHomogPlaneStressHomogenizer(c0,c1,dir,m1,frac);
            obj.lamTensor = lam.getPlaneStressHomogenizedTensor();
        end
        
        function createTolerance(obj)
           obj.tol = 1e-2; 
        end
    end
    
end

