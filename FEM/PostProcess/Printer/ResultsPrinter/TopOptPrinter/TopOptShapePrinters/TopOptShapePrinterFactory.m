classdef TopOptShapePrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(d,dT,shapeName)
            switch shapeName
                case 'ShFunc_NonSelfAdjoint_Compliance'
                    p = TopOptComplianceAndAdjointPrinter(d,dT);
                case 'ShFunc_Compliance'
                    p = TopOptCompliancePrinter(d,dT);
                case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
                    p = TopOptMicroPrinter(d,dT);
            end
        end
    end    
    
end