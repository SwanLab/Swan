classdef DualVariable < handle
    
    properties (Access = public)
       value 
       valueOld
    end
    
    properties (Access = private)
       nConstraints 
    end
    
   methods (Access = public, Static)
       
       function obj = DualVariable(cParams)
           obj.nConstraints = cParams.nConstraints;
           obj.value        = zeros(1,obj.nConstraints);           
       end
       
   end

    
end