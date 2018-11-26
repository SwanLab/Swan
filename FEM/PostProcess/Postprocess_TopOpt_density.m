classdef Postprocess_TopOpt_density < Postprocess_TopOpt
    
    properties (Access = private)
        FieldName = 'Density';
        ComponentName = 'Density'        
    end
    
    methods  (Access = protected)        
        
        function PrintResults(obj)
            obj.PrintScalar(obj.FieldName,obj.ComponentName,'Elastic Problem','Scalar','OnNodes','',obj.Field,obj.Iter);
        end
        
    end
    
    methods (Access = public)
        function getField(obj,results)
           obj.Field = results.density; 
        end
    end
end
