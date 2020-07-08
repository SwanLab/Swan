classdef SamplePointsCreatorForOptimalExponentComputerFactory < handle
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           
          switch cParams.type
              case 'FromMxMy'
                  obj = SamplePointsCreatorFromMxMyForOptimalExponentComputer(cParams);
              case 'FromFixedRho'
                  obj = SamplePointsCreatorFromFixedRhoForOptimalExponentComputer(cParams);
              case 'FromFixedRhoAndTxi'
                  obj = SamplePointsCreatorFromFixedRhoAndTxiForOptimalExponentComputer(cParams);
          end
           
       end
       
   end
   
end