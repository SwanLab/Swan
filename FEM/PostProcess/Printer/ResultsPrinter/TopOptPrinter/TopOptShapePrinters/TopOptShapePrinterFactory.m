classdef TopOptShapePrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(d,shapeName)
            switch shapeName
                case 'ShFunc_NonSelfAdjoint_Compliance'
                    p = TopOptComplianceAndAdjointPrinter(d);
                case 'ShFunc_Compliance'
                    p = TopOptCompliancePrinter(d);
                case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
                    p = TopOptMicroPrinter(d);
            end
        end
    end    
    
end