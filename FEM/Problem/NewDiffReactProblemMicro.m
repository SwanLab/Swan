classdef NewDiffReactProblemMicro < NewDiffReactProblem

    methods (Access = protected)
        
        function setScale(obj)
            obj.problemData.scale = 'MICRO';
        end
        
    end
    
end
