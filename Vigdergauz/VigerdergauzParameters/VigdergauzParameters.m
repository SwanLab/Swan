classdef VigdergauzParameters < handle 
    
   methods (Access = public, Static)
       
       function obj = create(cParams)
           f = VigdergauzParametersFactory();
           obj = f.create(cParams);
       end
       
   end
    
end