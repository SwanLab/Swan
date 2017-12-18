classdef Postprocess_TopOpt_levelSet < Postprocess_TopOpt_density
    
    
    properties
        levelSet_name = 'LevelSet';
        levelSet_name_component = 'LS';
    end
    
    methods (Access = public)
        
         function obj = Postprocess_TopOpt_levelSet()
         
         end
                    
        
    end
    
    
    methods (Access = protected)
        
        function Print_design_variable(obj,design_variable)
            obj.Print_LevelSet(design_variable);
        end
        
        function Print_LevelSet(obj,design_variable)
            obj.PrintScalar(obj.levelSet_name,obj.levelSet_name_component,'Elastic Problem','Scalar','OnNodes','',design_variable);
        end
        
    end
        
    
    
    
end