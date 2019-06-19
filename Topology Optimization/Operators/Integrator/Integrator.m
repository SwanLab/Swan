classdef Integrator < handle
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
    end
        
    methods (Access = protected, Abstract)        
    end
    
    methods (Access = public)

        
    end
    
    methods (Static, Access = public)
        
        function obj = create(cParams)
            obj = IntegratorFactory.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)           
           obj.mesh = cParams.mesh;
        end        
        
    end
    
    methods (Static, Access = protected)
        
        function quadrature = computeQuadrature(geometryType)
            quadrature = Quadrature.set(geometryType);
            quadrature.computeQuadrature('LINEAR');
        end        
        
    end
    
    methods (Access = private)
        
  
        
    end
    
    methods (Static, Access = private)
        
     
        
    end
    
end

