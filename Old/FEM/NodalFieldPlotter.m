classdef NodalFieldPlotter < handle
    
    properties (Access = private)
        mesh
        field
    end
    
    methods (Access = public)
        
        function obj = NodalFieldPlotter(cParams)
            obj.init(cParams)
        end
        
        function plot(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            z = obj.field;
            %figure()
            trisurf(obj.mesh.connec,x,y,z)
            view(0,90)
            colorbar
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.field = cParams.field;
        end
        
    end
    
end