classdef testHorizontalLaminate < TestSequentialLaminateTestedWithNumerics

    methods (Access = public)

        function obj = testHorizontalLaminate()
            obj.compute();
        end

    end
    
    methods (Access = protected)

        function loadLaminateDirection(obj)
            d = [0 1 0];
            dir = Vector3D;
            dir.setValue(d);
            dir.normalize()
            obj.LaminateDirection = dir;
        end

        function loadFiberDirection(obj)
            d = [1 0 0];
            dir = Vector3D;
            dir.setValue(d);
            dir.normalize()
            obj.FiberDirection = dir;
        end

    end

end