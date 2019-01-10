classdef OptimizerPrinter < handle
    
    properties (Access = private)
        postproc
        dataBase
        fields
        printMode
        mesh
    end
    
    methods (Access = public)
        
        function obj = OptimizerPrinter(mesh,optimizer,fileName,printMode,cost,constraint)
            obj.createPostProcessDataBaseStructre(mesh,optimizer,fileName,printMode,cost,constraint);
            obj.postproc = Postprocess('TopOptProblem',obj.dataBase);
            obj.mesh = mesh;
            obj.printMode = printMode;
        end
        
        function print(obj,x,iter,cost,constraint)
            obj.createTopOptFields(x,cost,constraint);
            obj.postproc.print(iter,obj.fields);
        end
        
    end
    
    methods (Access = private)
        
        function createPostProcessDataBaseStructre(obj,mesh,optimizer,fileName,printMode,cost,constraint)
            hasGaussData = false;
            quad = [];
            shapeFuncs = {cost.ShapeFuncs{:},constraint.ShapeFuncs{:}};
            nshapes = obtainNumberOfShapesToPrint(obj,shapeFuncs);
            shapesNames = cell(nshapes,1);
            iprint = 0;
            for ishape = 1:nshapes
                shapeName = class(shapeFuncs{ishape});
                if obj.hasToPrintShape(printMode,shapeName)
                    iprint = iprint + 1;
                    [hasGaussData,quad] = obj.hasGaussData(cost,constraint);
                    shapesNames{iprint} = shapeName;
                end
            end
            
            if isequal(printMode,'ElementalDensity') ...
                    || isequal(printMode,'DesignAndElementalDensity') ...
                    || isequal(printMode,'DesignElementalDensityAndShape')
                hasGaussData = true;
                phyPr = cost.ShapeFuncs{1}.getPhysicalProblem();
                quad = phyPr.element.quadrature;
            end
            
            d.quad = quad;
            d.iter = 0;
            d.mesh = mesh;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator.create(hasGaussData,d);
            obj.dataBase = ps.getValue();
            obj.dataBase.ShapeNames = shapesNames;
            obj.dataBase.optimizer  = optimizer;
            obj.dataBase.printMode  = printMode;
            obj.dataBase.hasGaussData = hasGaussData;
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
        
        function itHas = hasToPrintShape(obj,printMode,shape)
            itHas = false;
            switch  printMode
                case {'DesignAndShapes','DesignElementalDensityAndShape'}
                    switch shape
                        case {'ShFunc_NonSelfAdjoint_Compliance',...
                                'ShFunc_Compliance', ...
                                'ShFunc_Chomog_alphabeta',...
                                'ShFunc_Chomog_fraction'}
                            itHas = true;
                    end
                    
                    
            end
        end
        
        function [hasGaussData,quad] = hasGaussData(obj,cost,constraint)
            quad = [];
            hasGaussData = false;
            shapeFuncs = {cost.ShapeFuncs{:},constraint.ShapeFuncs{:}};
            shWithGaussData = {'ShFunc_NonSelfAdjoint_Compliance',...
                'ShFunc_Compliance',...
                'ShFunc_Chomog_alphabeta'};
            for ishape = 1:numel(shapeFuncs)
                shapeFun = shapeFuncs{ishape};
                isShapeWithGaussData = any(strcmp(shWithGaussData,class(shapeFun)));
                if isShapeWithGaussData
                    hasGaussData = true;
                    phyPr = shapeFun.getPhysicalProblem();
                    quad = phyPr.element.quadrature;
                    return
                end
            end
        end
        
        function createTopOptFields(obj,x,cost,constraint)
            
            if isequal(obj.printMode,'ElementalDensity') ...
                    || isequal(obj.printMode,'DesignAndElementalDensity')
                
                phi0 = obj.computePhiP0(x,cost);
                dens0  = obj.computeDensityP0(phi0);
                obj.fields.dens = dens0;
                
            elseif isequal(obj.printMode,'DesignElementalDensityAndShape')
                phi0 = obj.computePhiP0(x,cost);
                dens0  = obj.computeDensityP0(phi0);
                obj.fields.dens = dens0;
                
                shapeFuncs = {cost.ShapeFuncs{:},constraint.ShapeFuncs{:}};
                iprint = 0;
                for ishape = 1:numel(shapeFuncs)
                    shapeFun = shapeFuncs{ishape};
                    [f,hasToPrint] = obj.obtainFieldToPrint(shapeFun);
                    if hasToPrint
                        iprint = iprint + 1;
                        obj.fields.shVar{iprint} = f;
                    end
                end
                
            else
                
                shapeFuncs = {cost.ShapeFuncs{:},constraint.ShapeFuncs{:}};
                iprint = 0;
                for ishape = 1:numel(shapeFuncs)
                    shapeFun = shapeFuncs{ishape};
                    [f,hasToPrint] = obj.obtainFieldToPrint(shapeFun);
                    if hasToPrint
                        iprint = iprint + 1;
                        obj.fields.shVar{iprint} = f;
                    end
                end
                
            end
            obj.fields.designVariable = x;
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
        
        function phiP0 = computePhiP0(obj,x,cost)
            conec = obj.mesh.connec;
            phyPr = cost.ShapeFuncs{1}.getPhysicalProblem();
            shape = phyPr.element.interpolation_u.shape;
            quadr = phyPr.element.quadrature;
            ngaus = quadr.ngaus;
            nelem = size(conec,1);
            nnode = size(shape,1);
            
            phiP0 = zeros(ngaus,nelem);
            phi   = x;
            
            for igaus = 1:ngaus
                for inode = 1:nnode
                    nodes = conec(:,inode);
                    phiN(1,:) = phi(nodes);
                    phiP0(igaus,:) = phiP0(igaus,:) + shape(inode,igaus)*phiN;
                end
            end
            
        end
        
        function dens = computeDensityP0(obj,phi)
            dens = 1 - heaviside(phi);
        end
        
    end
    
end