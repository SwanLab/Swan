classdef Postprocess_TopOpt_levelSet < Postprocess_TopOpt
    properties
        levelSet_name = 'LevelSet';
        levelSet_name_component = 'LS';
    end
    
    methods (Access = public)        
        

        function PrintResults(obj)
            obj.PrintScalar(obj.levelSet_name,obj.levelSet_name_component,'Elastic Problem','Scalar','OnNodes','',obj.Field,obj.Iter);
        end
    end
    
    methods (Access = public)
        function getField(obj,results)
           obj.Field = results.design_variable; 
        end
    end
    
    
end