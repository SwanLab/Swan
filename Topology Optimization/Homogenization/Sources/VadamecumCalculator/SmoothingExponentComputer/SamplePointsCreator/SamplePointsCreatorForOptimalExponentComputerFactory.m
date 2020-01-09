classdef SamplePointsCreatorForOptimalExponentComputerFactory < handle
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           
          switch cParams.type
              case 'FromMxMy'
                  obj = SamplePointsCreatorFromMxMyForOptimalExponentComputer();
              case 'FromFixedRho'
                  obj = SamplePointsCreatorFromFixedRhoForOptimalExponentComputer(cParams);
          end
           
       end
       
   end
   
end