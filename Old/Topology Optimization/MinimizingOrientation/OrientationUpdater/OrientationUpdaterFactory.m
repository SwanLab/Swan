classdef OrientationUpdaterFactory < handle
        
   methods (Access = public, Static)
       
       function obj = create(cParams)
           
          switch cParams.type
              case 'MinimumEigenValue'
                  obj = OrientationUpdater_MinValue();
              case 'MaximumEigenValue'
                  obj = OrientationUpdater_MaxValue();
              case 'MinimumAbsEigenValue'
                  obj = OrientationUpdater_MinAbsValue();
              case 'MaximumAbsEigenValue'
                  obj = OrientationUpdater_MaxAbsValue();
          end
           
       end
       
   end
    
end