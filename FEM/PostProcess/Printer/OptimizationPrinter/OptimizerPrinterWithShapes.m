classdef OptimizerPrinterWithShapes <  OptimizerPrinter
    
    properties (Access = private)
       shapeNames
       shapeFuncs       
       quad
    end
    
    methods (Access = protected)
        
        function obtainHasGaussDataAndQuad(obj)
            obj.shapeFuncs = {obj.cost.ShapeFuncs{:},obj.constraint.ShapeFuncs{:}};            
            h = false;
            q = [];
            s = [];
            if obj.thereAreShapesToPrint(obj.shapeFuncs)
                nshapes = obj.obtainNumberOfShapesToPrint(obj.shapeFuncs);
                [sN,sF] = obj.obtainPrintableShapes(obj.shapeFuncs,nshapes);
                s = sN;
                [h,q] = obj.obtainQuadHasGaussData(sF{1});
            end
            obj.shapeNames = s;
            obj.hasGaussData = h;
            obj.quad = q;
            
            if obj.hasGaussData
                obj.dI.quad = obj.quad;
            end
        end
        
        function createDataBaseForPostProcess(obj,ps,optimizer,printMode)
            obj.dataBase = ps.getValue();
            obj.dataBase.hasGaussData = obj.hasGaussData;            
            obj.dT.optimizer = optimizer;
            obj.dT.printMode = printMode;
            obj.dT.ShapeNames = obj.shapeNames;                   
        end        
        
    end
    
    
    
    methods (Access = private)
        
        function itHas = thereAreShapesToPrint(obj,shapeFuncs)
            itHas = false;
            for ishape = 1:numel(shapeFuncs)
                shapeFunc = shapeFuncs{ishape};
                if obj.isShapePrintable(shapeFunc)
                    itHas = true;
                    return
                end
            end
        end
        
        function n = obtainNumberOfShapesToPrint(obj,shapeFuncs)
            n = 0;
            for ishape = 1:numel(shapeFuncs)
                shapeName = class(shapeFuncs{ishape});
                switch shapeName
                    case {'ShFunc_NonSelfAdjoint_Compliance',...
                            'ShFunc_Compliance', ...
                            'ShFunc_Chomog_alphabeta',...
                            'ShFunc_Chomog_fraction'}
                        n = n +1;
                end
            end
        end
        
        function [sN,sF] = obtainPrintableShapes(obj,shapeFuncs,nshapes)
            shapesNamesToPrint = cell(nshapes,1);
            shapeFuncsToPrint = cell(nshapes,1);
            iprint = 0;
            for ishape = 1:nshapes
                shapeFunc = shapeFuncs{ishape};
                shapeName = class(shapeFunc);
                if obj.isShapePrintable(shapeFunc)
                    iprint = iprint + 1;
                    shapesNamesToPrint{iprint} = shapeName;
                    shapeFuncsToPrint{iprint}  = shapeFunc;
                end
            end
            sN = shapesNamesToPrint;
            sF = shapeFuncsToPrint;
        end
        
    end
    
    
end