classdef Geometry < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = GeometryFactory();
            obj = f.create(cParams);
        end               
        
    end
    
end