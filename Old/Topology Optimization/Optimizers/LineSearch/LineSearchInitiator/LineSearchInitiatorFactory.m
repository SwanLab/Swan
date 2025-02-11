classdef LineSearchInitiatorFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'STANDARD'
                    switch cParams.optimizerType
                        case 'PROJECTED GRADIENT'
                            obj = Normalized_LineSearchInitiator(cParams);
                        case {'SLERP','HAMILTON-JACOBI'}
                            cParams.value = 1;
                            obj = Constant_LineSearchInitiator(cParams);
                    end
                case 'INCREASING LAST STEP'
                    switch cParams.optimizerType                        
                        case {'PROJECTED GRADIENT','HAMILTON-JACOBI'}
                            cParams.maxValue = +Inf;
                        case {'SLERP'}
                            cParams.maxValue = 1;
                    end
                    obj = IncreaseLast_LineSearchInitiator(cParams);
            end
            
        end
        
    end
    
end

