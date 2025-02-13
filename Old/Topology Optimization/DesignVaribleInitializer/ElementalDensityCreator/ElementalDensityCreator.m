classdef ElementalDensityCreator < handle
    
    
   properties (Access = protected)
       density
   end
   
   methods (Access = public, Static)
       
       function eC = create(type,cParams)
          f = ElementalDensityCreatorFactory();
          eC = f.create(type);
          eC.createDensity(cParams);
       end
   end
   
   methods (Access = public)
       
       function dC = getDensityCreator(obj)
           dC = obj.densityCreator;
       end
       
       function d = getDensity(obj)
           d = obj.density;
       end
       
       function f = getFieldsToPrint(obj)
           f{1} = obj.density;
       end
       
   end
    
   methods (Access = protected, Abstract)
       createDensity(obj)
   end
    
end