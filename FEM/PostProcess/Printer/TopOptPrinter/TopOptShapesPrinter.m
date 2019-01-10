classdef TopOptShapesPrinter < TopOptPrinter
    
    properties (Access = private)
        ndim
    end
    
    methods (Access = public)
        
        function obj = TopOptShapesPrinter(ndim)
            obj.ndim = ndim;
        end
        
        function create(obj,d)
            shapeFunNames = d.ShapeNames;
            simStr = obj.simulationStr;
            for ishape = 1:numel(shapeFunNames)
                switch shapeFunNames{ishape}
                    case 'ShFunc_NonSelfAdjoint_Compliance'
                        p = TopOptComplianceAndAdjointPrinter(simStr);
                    case 'ShFunc_Compliance'
                        p = TopOptCompliancePrinter(simStr);
                    case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
                        p = TopOptMicroPrinter(simStr,obj.ndim);
                end
                shPrint = p.createShapePrinters(d);
                obj.printers{ishape} = shPrint;
            end
        end
        
        function print(obj,istep,fields)
            is = istep;
            for ishape = 1:numel(fields.shVar)
                shVar     = fields.shVar{ishape};
                shPrinter = obj.printers{ishape};
                for ivar = 1:numel(shVar)
                    var = shVar{ivar};
                    printer = shPrinter{ivar};
                    printer.printOnlyResults(is,var);
                end
            end
        end
        
    end
    
end