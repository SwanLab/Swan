classdef testMaterialInsertingLpBall < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testMaterialInsertingLpBall';
    end
    
    properties (Access = private)
        vadSAcomparator
        plotter
    end
    
    
    methods (Access = public)
        
        function obj = testMaterialInsertingLpBall()
            obj.computeVademecumAndSimpAllValues();
            obj.plotVademecumAndSimpAllValues();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function computeVademecumAndSimpAllValues(obj)
            s.fileName =  'VademecumSmoothCorner';
            v = VademecumSimpAllComparator(s);
            v.calculate();
            obj.vadSAcomparator = v;
        end
        
        function plotVademecumAndSimpAllValues(obj)
            s.vadSAcomparator = obj.vadSAcomparator;
            obj.plotter = VademecumComparatorPlotter(s);
            obj.plotter.plot();
        end
        
      
        
    end    
    
end