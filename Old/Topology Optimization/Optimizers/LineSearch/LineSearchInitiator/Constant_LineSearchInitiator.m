classdef Constant_LineSearchInitiator < LineSearchInitiator
    
    properties (Access = private)
        value
    end

    methods (Access = public)
        
        function obj = Constant_LineSearchInitiator(cParams)
            obj.init(cParams);
            obj.value = cParams.value;
        end
        
        function initStep = compute(obj,lastStep)            
            initStep = obj.value;
        end       
        
    end     
  
  
end

