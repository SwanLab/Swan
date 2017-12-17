classdef Postprocess_TopOpt < Postprocess_PhysicalProblem
    
    properties  

    end
    
    
    methods (Access = public)
        function obj = Postprocess_TopOpt()
            

            
        end
    end
        
    methods (Static)    
        
        function obj = Create(optimizer)
            switch optimizer
                case 'SLERP'
                    obj = Postprocess_TopOpt_levelSet();
                case 'PROJECTED GRADIENT'
                    obj = Postprocess_TopOpt_density();
                case 'MMA'
                    obj = Postprocess_TopOpt_density();
            end            
            
        end
        
    end
end
