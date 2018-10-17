classdef PostprocessDensityInGaussPoints < Postprocess_TopOpt
    
    
    properties (Access = private)
        FieldName = 'RegularizedDensity';
        ComponentName = 'RegularizedDensity'
    end
    
    methods (Access = public)
        
        function obj = PostprocessDensityInGaussPoints(quadrature)
            obj.SetQuadratureInfo(quadrature)
        end
        
        function PrintResults(obj)
            obj.ngaus = size(obj.Field,2);
            obj.PrintGaussPointsHeader()
            obj.PrintScalar(obj.FieldName,obj.ComponentName,'Elastic Problem','Scalar','OnGaussPoints',obj.gauss_points_name,obj.Field,obj.Iter);
        end
       
    end
    
    methods (Access = public)
        function getField(obj,results)
            obj.Field = results.density;
        end
    end
    
    methods (Access = private)
        function SetQuadratureInfo(obj,quadrature)
           obj.ngaus = quadrature.ngaus;
           obj.posgp = quadrature.posgp';        
        end
    end
    
end

