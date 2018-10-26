classdef testHorizontalTensorRotatedVsRank2 < ...
        testHorizontalTensorRotatedVsSequentialLaminate
    
    properties (Access = protected)
        tol
    end
    
    methods (Access = public)
        
        function obj = testHorizontalTensorRotatedVsRank2()
            obj.computeTest()
        end
    end
    
    methods (Access = protected)
        
        function computeLaminateDirectly(obj)
            c0    = obj.C0;
            c0.computeTensorVoigt()
            c0.computeTensorVoigtInPlaneStress()
            c1    = obj.C1;
            c1.computeTensorVoigt()
            c1.computeTensorVoigtInPlaneStress()                                    
            dir(1,:) = obj.lamDir;
            dir(2,:) = obj.lamDir;
            mi(1)    = obj.lamPar;
            mi(2)    = 0;
            frac  = obj.theta;
            rank2 = RankTwoLaminateHomogenizer(c0,c1,dir,mi,frac);
            obj.lamTensor = rank2.Ch;
        end
        
        function createTolerance(obj)
           obj.tol = 1e-2; 
        end
    end
    
end

