classdef LineSearch_Backtracking < LineSearch
    
    properties (Access = private)
       rate 
    end
    
    methods (Access = public)
        
        function obj = LineSearch_Backtracking(cParams)
            obj.init(cParams);
            obj.rate = cParams.rate;
        end
        
        function update(obj)
            obj.value = obj.value*obj.rate;
            obj.nTrials = obj.nTrials + 1;
        end        
        
    end   
    
end