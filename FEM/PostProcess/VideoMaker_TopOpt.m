classdef VideoMaker_TopOpt < VideoMaker_Physical_Problem
    
    
     properties 
         
         
     end
       

    methods (Access = public)
        function obj = VideoMaker_TopOpt()
            

        end
    end
        
   methods (Static)
       
       function obj = Create(optimizer,pdim)
           switch pdim
               case '2D'
                   switch optimizer
                       case 'SLERP'
                           obj = VideoMaker_TopOpt_levelSet();
                       case 'PROJECTED GRADIENT'
                           obj = VideoMaker_TopOpt_density();
                       case 'MMA'
                           obj = VideoMaker_TopOpt_density();
                       case 'IPOPOT'
                           obj = VideoMaker_TopOpt_density();
                   end
               case '3D'
                    switch optimizer
                       case 'SLERP'
                           obj = VideoMaker_TopOpt_levelSet3D();
                       case 'PROJECTED GRADIENT'
                           obj = VideoMaker_TopOpt_density3D();
                       case 'MMA'
                           obj = VideoMaker_TopOpt_density3D();
                       case 'IPOPOT'
                           obj = VideoMaker_TopOpt_density3D();
                    end
           end            
        end
    end

end