classdef OptimizerPrinterWithShapes <  OptimizerPrinter
    
    properties (Access = private)
        shapeNames
        shapeFuncs
        quad
    end
    
    methods (Access = protected)
        
        function createDtDataBase(obj,optimizer,printMode)
            obj.createDtDataBase@OptimizerPrinter(optimizer,printMode)
            obj.dT.ShapeNames = obj.shapeNames;
        end
        
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
        end
        
        function createTopOptFields(obj,x,cost,constraint)
            shFuncs = {cost.ShapeFuncs{:},constraint.ShapeFuncs{:}};
            iprint = 0;
            for ishape = 1:numel(shFuncs)
                shapeFun = shFuncs{ishape};
                [f,hasToPrint] = obj.obtainFieldToPrint(shapeFun);
                if hasToPrint
                    iprint = iprint + 1;
                    obj.fields.shVar{iprint} = f;
                end
            end
        end
        
        function createDataInputForCreateDataBase(obj,mesh,fileName)
            obj.createDataInputForCreateDataBase@OptimizerPrinter(mesh,fileName);
            if obj.hasGaussData
               obj.dI.quad = obj.quad;
            end
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
        
        function itIs = isShapePrintable(obj,shapeFunc)
            shapeName = class(shapeFunc);
            printingShapes = {'ShFunc_NonSelfAdjoint_Compliance',...
                'ShFunc_Compliance', ...
                'ShFunc_Chomog_alphabeta',...
                'ShFunc_Chomog_fraction'};
            itIs = any(strcmp(printingShapes,shapeName));
        end
        
        function [f,hasToPrint] = obtainFieldToPrint(obj,shapeFun)
            hasToPrint = false;
            f = [];
            shapeFunName = class(shapeFun);
            switch shapeFunName
                case 'ShFunc_NonSelfAdjoint_Compliance'
                    hasToPrint = true;
                    phyPr = shapeFun.getPhysicalProblem();
                    f{1} = phyPr.variables;
                    adjPr = shapeFun.getAdjointProblem();
                    f{2} = adjPr.variables;
                case 'ShFunc_Compliance'
                    hasToPrint = true;
                    phyPr = shapeFun.getPhysicalProblem();
                    f{1} = phyPr.variables;
                case {'ShFunc_Chomog_alphabeta','ShFunc_Chomog_fraction'}
                    hasToPrint = true;
                    phyPr = shapeFun.getPhysicalProblem();
                    var = phyPr.variables.var2print;
                    f = cell(numel(var),1);
                    for it = 1:numel(var)
                        f{it} = var{it};
                    end
            end
        end        
        
        function [h,q] = obtainQuadHasGaussData(obj,sh)
            h = true;
            phyPr = sh.getPhysicalProblem();
            q = phyPr.element.quadrature;
        end
    end
    
    
    
end