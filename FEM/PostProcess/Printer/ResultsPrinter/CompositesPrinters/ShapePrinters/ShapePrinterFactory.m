classdef ShapePrinterFactory < handle
    
    
    methods (Access = public, Static)
        
        function p = create(d,shapeName)
            switch shapeName
                case 'ShFunc_NonSelfAdjoint_Compliance'
                    p = ComplianceAndAdjointPrinter(d);
                case 'ShFunc_Compliance'
                    p = CompliancePrinter(d);
                case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
                    p = HomogenizedTensorPrinter(d);
            end
        end
    end    
    
end