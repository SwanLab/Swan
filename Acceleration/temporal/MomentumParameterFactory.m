classdef MomentumParameterFactory < handle
   
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'CONSTANT'
                    obj = MomentumParameterConstant(cParams);
                case 'CONVEX'
                    obj = MomentumParameterConvexCase(cParams);
                case 'STRICTLY CONVEX'
            
            end
         
        end
        
    end

end