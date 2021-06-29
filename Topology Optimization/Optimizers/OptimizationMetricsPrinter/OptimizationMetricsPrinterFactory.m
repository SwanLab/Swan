classdef OptimizationMetricsPrinterFactory < handle
    
    methods (Access = public, Static)
        
        function printer = create(cParams)
            if cParams.shallPrint
                switch cParams.optimizer.type
                    case {'AlternatingPrimalDual','DualNestedInPrimal'}
                        printer = OptimizationMetricsPrinter_AugLag(cParams);
                    case 'MMA'
                        printer = OptimizationMetricsPrinter_MMA(cParams);
                    case 'IPOPT'
                        printer = OptimizationMetricsPrinter_IPOPT(cParams);
                    otherwise
                        error('Invalid optimizer type.')
                end
            else
                printer = OptimizationMetricsPrinter_Null(cParams);
            end
        end
        
    end
    
end