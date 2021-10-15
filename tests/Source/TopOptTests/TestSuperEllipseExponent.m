classdef TestSuperEllipseExponent < handle
    
    properties (Access = private)
        dataFile = 'test_superEllipseExponent'
        m1
        m2
        qD
        q
    end
    
    methods (Access = public)
        
        function obj = TestSuperEllipseExponent()
            obj.loadData();
            obj.computeExponent();
        end

        function error = computeError(obj)
            error = norm(obj.qD - obj.q(:))/norm(obj.qD);
        end

    end
    
    methods (Access = private)
        
        function loadData(obj)
            d = load(obj.dataFile);
            obj.m1 = d.m1;
            obj.m2 = d.m2;
            obj.qD = d.q;
        end
        
        function computeExponent(obj)
            s.m1 = obj.m1;
            s.m2 = obj.m2;
            o = SmoothingExponentComputerOptimal(s);
            obj.q = o.compute();
        end

    end

end