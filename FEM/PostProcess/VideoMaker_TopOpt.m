classdef VideoMaker_TopOpt < VideoMaker_Physical_Problem
    
    
     properties 
         
         
     end
       

    methods (Access = public)
        function obj = VideoMaker_TopOpt()
            

        end
    end
        
   methods (Static)
        
        function obj = Create(optimizer)
            switch optimizer
                case 'SLERP'
                    obj = VideoMaker_TopOpt_levelSet();
                case 'PROJECTED GRADIENT'
                    obj = VideoMaker_TopOpt_density();
                case 'MMA'
                    obj = VideoMaker_TopOpt_density();
            end
            
        end
    end

end