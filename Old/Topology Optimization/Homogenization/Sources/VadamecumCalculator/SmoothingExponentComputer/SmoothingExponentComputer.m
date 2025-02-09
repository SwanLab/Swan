classdef SmoothingExponentComputer < handle

    properties (Access = protected)
        value
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = SmoothingExponentComputerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
                       
        function q = compute(obj)
           obj.computeExponent(); 
           q = obj.value;
        end 
        
    end
    
    methods (Access = protected, Abstract)
       computeExponent(obj); 
    end

end