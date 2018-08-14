classdef Postprocess_TopOpt_levelSet < Postprocess_TopOpt_density
    properties
        levelSet_name = 'LevelSet';
        levelSet_name_component = 'LS';
    end
    
    methods (Access = public)        
        function Print_design_variable(obj,design_variable)
            obj.Print_LevelSet(design_variable);
        end
        
        function Print_LevelSet(obj,results)
            obj.PrintScalar(obj.levelSet_name,obj.levelSet_name_component,'Elastic Problem','Scalar','OnNodes','',results.design_variable,results.iter);
        end
    end
end