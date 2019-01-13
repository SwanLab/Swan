classdef TopOptShapesPrinter < TopOptPrinter
            
    methods (Access = public)
       
        function obj = TopOptShapesPrinter(d,dT,s)         
            shapeNames = s;
            nShapes = numel(shapeNames);
            obj.printers = cell(nShapes,1);
            for ishape = 1:nShapes
                shapeName = shapeNames{ishape};
                p = TopOptShapePrinter.create(d,dT,shapeName);   
                obj.printers{ishape} = p.getPrinters();
            end
        end

        function print(obj,istep,fields)
            is = istep;
            for ishape = 1:numel(fields.shVar)
                shVar     = fields.shVar{ishape};
                shPrinter = obj.printers{ishape};
                for ivar = 1:numel(shVar)
                    var = shVar{ivar};
                    p = shPrinter{ivar};
                    p.printOnlyResults(is,var);
                end
            end
        end
        
    end
    
end