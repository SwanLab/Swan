classdef PlaneStressIndex < handle
    
    properties (Access = private)
        InPlane
        OutPlane
    end
    
    methods (Access = public)
        
        function obj = PlaneStressIndex()
            obj.createInPlaneIndex()
            obj.createOutPlaneIndex()
        end
        
        function Index = getInPlaneIndex(obj)
            Index = obj.InPlane;
        end
        
        function Index = getOutPlaneIndex(obj)
            Index = obj.OutPlane;
        end
    end
    
    methods (Access = private)
        
        function createInPlaneIndex(obj)
            obj.InPlane = [1 2 6];
        end
        
        function createOutPlaneIndex(obj)
            obj.OutPlane = [3 4 5];
        end
    end
    
end

