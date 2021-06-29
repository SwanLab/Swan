classdef Geometry < handle
    
    properties (GetAccess = public, SetAccess = protected)
        dvolu
    end    
    
    properties (SetAccess = private, GetAccess = protected)
        quadrature
        interpolationVariable
        mesh
    end
    
    properties (Access = protected)
        coordElem            
        permutation
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = GeometryFactory();
            obj = f.create(cParams);
        end               
        
    end
    
    methods (Access = protected)
        
        function initGeometry(obj,interpV,quad)
            obj.interpolationVariable = interpV;
            obj.quadrature = quad;
            obj.computeShapeFunctions();
        end        
       
       function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.coordElem = permute(obj.mesh.coordElem,obj.permutation);                         
       end        
       
    end
    
    methods (Access = private)
        
        function computeShapeFunctions(obj)
            xpg = obj.quadrature.posgp;
            obj.interpolationVariable.computeShapeDeriv(xpg)
            obj.mesh.interpolation.computeShapeDeriv(xpg);
        end        
       
    end
    
end