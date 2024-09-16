classdef Optimizer < handle

    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        monitoring
        maxIter
        nIter = 0
        tolerance
        dualUpdater
        primalUpdater
        constraintCase
        postProcess
    end

    properties (Access = public)
        simulationPrinter
    end

    properties (Access = private)
        outFilename % !!!
        outFolder
    end


    methods (Access = public, Static)

        function obj = create(cParams)
            f   = OptimizerFactory();
            obj = f.create(cParams);
        end

    end

    methods (Access = protected)

        function initOptimizer(obj,cParams)
            obj.nIter          = 0;
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.maxIter        = cParams.maxIter;
            obj.tolerance      = cParams.tolerance;
            obj.constraintCase = cParams.constraintCase;
            %obj.createPostProcess(cParams.postProcessSettings);
        end

        function createPrimalUpdater(obj,cParams)
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(cParams);
        end

        function createDualUpdater(obj,cParams)
            cParams.type    = obj.type;
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);
        end

        function isAcceptable = checkConstraint(obj)
            for i = 1:length(obj.constraint.value)
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        fine = obj.checkEqualityConstraint(i);
                    case {'INEQUALITY'}
                        fine = obj.checkInequalityConstraint(i);
                end
                if fine
                    isAcceptable = true;
                else
                    isAcceptable = false;
                    break;
                end
            end
        end

        function printOptimizerVariable(obj)
            if ~isempty(obj.postProcess)
                d.fields  = obj.designVariable.getVariablesToPlot();
                d.cost = obj.cost;
                d.constraint = obj.constraint;
                %                 obj.postProcess.print(obj.nIter,d);
                [desFun, desName] = obj.designVariable.getFunsToPlot();
                fun  = desFun;
                name = desName;
                for iShp = 1:numel(obj.cost.shapeFunctions)
                    [shpFun, shpName] = obj.cost.shapeFunctions{iShp}.getFunsToPlot();
                    fun  = [fun, shpFun];
                    name = [name, shpName];
                end
                file = [obj.outFolder,'/',obj.outFilename, '_', num2str(obj.nIter)];

                zz.mesh     = obj.designVariable.mesh;
                zz.filename = file;
                zz.fun      = fun;
                zz.funNames = name;
                pp = FunctionPrinter_Paraview(zz);
                pp.print();
                obj.simulationPrinter.appendStep(file);
            end
            obj.obtainGIF();
            if ismethod(obj.designVariable,'plot')
                obj.designVariable.plot();
            end
        end

        function obtainGIF(obj)
            %set(0,'DefaultFigureVisible','off');






            gifName = 'NullSLERPResults/TopOpt/MultiLoadBridge/NoOscillations/gJ10_9Loads';







            deltaTime = 0.01;
            m = obj.designVariable.fun.mesh;
            xmin = min(m.coord(:,1));
            xmax = max(m.coord(:,1));
            ymin = min(m.coord(:,2));
            ymax = max(m.coord(:,2));

            f = obj.designVariable.fun.fValues;
            switch obj.designVariable.type
                case 'LevelSet'
                    uMesh = obj.designVariable.getUnfittedMesh();
                    uMesh.compute(f);
                    gifFig = figure;
                    uMesh.plotStructureInColor('black');
                    hold on
                case 'Density'
                    p1.mesh    = m;
                    p1.fValues = f;
                    p1.order   = 'P1';
                    RhoNodal   = LagrangianFunction(p1);
                    q = Quadrature.set(m.type);
                    q.computeQuadrature('CONSTANT');
                    xV = q.posgp;
                    RhoElem = squeeze(RhoNodal.evaluate(xV));

                    figHandle = figure;
                    axis off
                    axis equal
                    axes = figHandle.Children;
                    patchHandle = patch(axes,'Faces',m.connec,'Vertices',m.coord,...
                        'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
                    set(axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
                    set(patchHandle,'FaceVertexAlphaData',RhoElem,'FaceAlpha','flat');
            end
            fig = gifFig;
            fig.CurrentAxes.XLim = [xmin xmax];
            fig.CurrentAxes.YLim = [ymin ymax];
            axis([xmin xmax ymin ymax])
            gifname = [gifName,'.gif'];
            set(gca, 'Visible', 'off')

            frame = getframe(fig);
            [A,map] = rgb2ind(frame.cdata,256);
            if obj.nIter == 0
                imwrite(A,map,gifname,"gif","LoopCount",0,"DelayTime",deltaTime);
            else
                imwrite(A,map,gifname,"gif","WriteMode","append","DelayTime",deltaTime);
            end
            close(gifFig);

            %set(0,'DefaultFigureVisible','on');
        end

    end

    methods (Access = private)
        function createPostProcess(obj,cParams)
            if cParams.shallPrint
                d = obj.createPostProcessDataBase(cParams);
                d.printMode = cParams.printMode;
                d.nDesignVariables = obj.designVariable.nVariables;
                obj.postProcess = Postprocess('TopOptProblem',d);
                s.filename = [obj.outFolder,'/',obj.outFilename, '_simulation'];
                obj.simulationPrinter = SimulationPrinter(s);
            end
        end

        function d = createPostProcessDataBase(obj,cParams)
            path = pwd;
            obj.outFolder   = fullfile(path,'Output',cParams.femFileName);
            obj.outFilename = cParams.femFileName;
            d.mesh    = obj.designVariable.mesh;
            d.outFileName = cParams.femFileName;
            d.ptype   = cParams.ptype;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.create();
            d.ndim       = obj.designVariable.mesh.ndim;
            d.pdim       = cParams.pdim;
            d.optimizer  = obj.type;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;
        end

        function c = checkInequalityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = g <= 0;
        end

        function c = checkEqualityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = abs(g) < obj.tolerance;
        end

    end

end