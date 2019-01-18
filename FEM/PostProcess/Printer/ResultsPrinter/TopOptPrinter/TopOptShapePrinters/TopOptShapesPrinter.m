classdef TopOptShapesPrinter < TopOptPrinter
    
    properties (Access = private)
        allShapes
        printableIndex
    end
    
    methods (Access = public)
        
        function obj = TopOptShapesPrinter(d)
            obj.obtainAllShapes(d);
            obj.obtainPrintableIndex();
            obj.createShapesPrinters(d);
        end
        
        function print(obj,istep)
            for iprinter = 1:numel(obj.printers)
                p = obj.printers{iprinter};
                p.print(istep)
            end  
        end
        
        function itHas = hasGaussData(obj)
            itHas = false;
            for iprinter = 1:numel(obj.printers)
                hasGaussDataPrinter = obj.printers{iprinter}.hasGaussData();
                itHas = hasGaussDataPrinter || itHas;
            end
        end
        
        function storeResultsInfo(obj,d)
            obj.obtainAllShapes(d);
            obj.obtainPrintableIndex();                        
            for iprinter = 1:numel(obj.printers)
                index = obj.printableIndex(iprinter);
                pS = obj.allShapes{index};
                p = obj.printers{iprinter};
                p.storeResultsInfo(pS);
            end  
        end
        
    end

    
    methods (Access = private)
        
        function obtainAllShapes(obj,d)
            obj.allShapes = {d.cost.ShapeFuncs{:},d.constraint.ShapeFuncs{:}};
        end
        
        function obtainPrintableIndex(obj)
            aShapes = obj.allShapes;
            nAllShapes = numel(aShapes);
            iprint = 0;
            for ishape = 1:nAllShapes
                shape = aShapes{ishape};
                if obj.isShapePrintable(shape)
                    iprint = iprint + 1;
                    obj.printableIndex(iprint) = ishape;
                end
            end
        end
        
        function createShapesPrinters(obj,d)
            pI = obj.printableIndex;
            aS = obj.allShapes;
            nShapes = numel(pI);
            obj.printers = cell(nShapes,1);
            factory = TopOptShapePrinterFactory();
            for ishape = 1:nShapes
                index = obj.printableIndex(ishape);
                shapeName = class(aS{index});
                p = factory.create(d,shapeName);
                obj.printers{ishape} = p;
            end
        end

        
%         function [f,hasToPrint] = obtainFieldToPrint(obj,shapeFun)
%             hasToPrint = false;
%             f = [];
%             shapeFunName = class(shapeFun);
%             switch shapeFunName
%                 case 'ShFunc_NonSelfAdjoint_Compliance'
%                     hasToPrint = true;
%                     phyPr = shapeFun.getPhysicalProblem();
%                     f{1} = phyPr.variables;
%                     adjPr = shapeFun.getAdjointProblem();
%                     f{2} = adjPr.variables;
%                 case 'ShFunc_Compliance'
%                     hasToPrint = true;
%                     phyPr = shapeFun.getPhysicalProblem();
%                     f{1} = phyPr.variables;
%                 case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
%                     hasToPrint = true;
%                     phyPr = shapeFun.getPhysicalProblem();
%                     var = phyPr.variables.var2print;
%                     f = cell(numel(var),1);
%                     for it = 1:numel(var)
%                         f{it} = var{it};
%                     end
%             end
%         end
%         
%         function itHas = thereAreShapesToPrint(obj,shapeFuncs)
%             itHas = false;
%             for ishape = 1:numel(shapeFuncs)
%                 shapeFunc = shapeFuncs{ishape};
%                 if obj.isShapePrintable(shapeFunc)
%                     itHas = true;
%                     return
%                 end
%             end
%         end
        
        
        function itIs = isShapePrintable(obj,shapeFunc)
            shapeName = class(shapeFunc);
            printingShapes = {'ShFunc_NonSelfAdjoint_Compliance',...
                'ShFunc_Compliance', ...
                'ShFunc_Chomog_alphabeta',...
                'ShFunc_Chomog_fraction'};
            itIs = any(strcmp(printingShapes,shapeName));
        end
        
    end
    
end