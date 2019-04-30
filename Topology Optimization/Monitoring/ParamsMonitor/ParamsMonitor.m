classdef ParamsMonitor < ParamsMonitor_Interface
    
    properties (Access = private)
        frame
        figures
        nFigs
        nRows
        nCols
        
        problemID
        costFuncValue
        nCostFuncs
        constraintValues
        nConstraints
        refreshInterval
        convergenceVars
        iteration
        hasFinished
        iRefresh
        iStep
        nStep
        
        namingManager
    end
    
    methods (Access = public)
        
        function obj = ParamsMonitor(cParams)
            obj.init(cParams);
            drawnow
        end
        
        function refresh(obj,it,hasFinished,iStep,nStep)
            obj.updateParams(it,hasFinished,iStep,nStep);
            obj.updateFigures();
            if obj.shallRefresh()
                obj.refreshFigures();
                obj.refreshFrameTitle();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.problemID = cParams.problemID;
            obj.refreshInterval = cParams.refreshInterval;
            obj.nCostFuncs = length(cParams.costFuncNames);
            obj.nConstraints = length(cParams.constraintFuncs);
            
            obj.costFuncValue    = cParams.cost;
            obj.constraintValues = cParams.constraint;
            obj.convergenceVars  = cParams.convergenceVars;
            obj.createNamingManager(cParams);
            obj.createFigures();
        end
        
        function createFigures(obj)
            obj.figures = {};
            obj.initCostFigures();
            obj.initConstraintFigures();
            obj.initConvergenceVarsFigures();
            obj.displayFigures();
        end
        
        function initCostFigures(obj)
            obj.initFigure('Cost');
            for i = 1:obj.nCostFuncs
                title = obj.namingManager.getCostFuncFigureTitle(i);
                obj.initFigure(title);
            end
        end
        
        function initConstraintFigures(obj)
            for i = 1:obj.nConstraints
                obj.initConstraintFigure(i);
                obj.initLambdaFigure(i);
            end
        end
        
        function initConstraintFigure(obj,i)
            title = obj.namingManager.getConstraintFigureTitle(i);
            obj.initFigure(title);
        end
        
        function initLambdaFigure(obj,i)
            title = obj.namingManager.getLambdaFigureTitle(i);
            obj.initFigure(title);
        end
        
        function initConvergenceVarsFigures(obj)
            for i = 1:obj.convergenceVars.nVar
                title = obj.namingManager.getConvVarFigureTitle(i);
                obj.initFigure(title);
            end
        end
        
        function obj = initFigure(obj,title)
            chartType = obj.getChartType(title);
            newFig = DisplayFactory.create(chartType,title);
            
            obj.appendFigure(newFig);
        end
        
        function appendFigure(obj,fig)
            obj.figures{end+1} = fig;
        end
        
        function displayFigures(obj)
            obj.createFrame();
            obj.computeDistribution();
            for i = 1:obj.nFigs
                obj.figures{i}.show(obj.nRows,obj.nCols,i,[0.06 0.03]);
            end
        end
        
        function createFrame(obj)
            obj.frame = figure;
            warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
            drawnow; set(get(obj.frame,'JavaFrame'),'Maximized',1);
        end
        
        function computeDistribution(obj)
            if obj.nFigs <= 4
                obj.nRows = 1;
            else
                obj.nRows = 2;
            end
            obj.nCols = round(obj.nFigs/obj.nRows);
        end
        
        function updateParams(obj,it,hasFinished,iStep,nStep)
            obj.iteration = it;
            obj.hasFinished = hasFinished;
            obj.iRefresh = 1;
            obj.iStep = iStep;
            obj.nStep = nStep;
        end
        
        function updateFigures(obj)
            obj.updateCost();
            obj.updateConstraints();
            obj.updateConvergenceVars();
        end
        
        function updateCost(obj)
            obj.updateFigure(obj.costFuncValue.value);
            for i = 1:obj.nCostFuncs
                obj.updateFigure(obj.costFuncValue.shapeFunctions{i}.value)
            end
        end
        
        function updateConstraints(obj)
            for i = 1:obj.nConstraints
                obj.updateFigure(obj.constraintValues.shapeFunctions{i}.value)
                if ~isempty(obj.constraintValues.lambda)
                    obj.updateFigure(obj.constraintValues.lambda(i))
                end
            end
        end
        
        function updateConvergenceVars(obj)
            for i = 1:obj.convergenceVars.nVar
                if ~isempty(obj.convergenceVars.value(i))
                    obj.updateFigure(obj.convergenceVars.value(i));
                end
            end
        end
        
        function updateFigure(obj,value)
            obj.figures{obj.iRefresh}.updateParams(obj.iteration,value,obj);
        end
        
        function refreshFigures(obj)
            for iFig = 1:obj.nFigs
                obj.figures{iFig}.refresh();
            end
        end
        
        function refreshFrameTitle(obj)
            title = obj.namingManager.getFrameTitle(obj.iStep,obj.nStep,obj.iteration);
            set(obj.frame,'NumberTitle','off','Name',title)
        end
        
        function itShall = shallRefresh(obj)
            itShall = (obj.isTimeToRefresh() || obj.hasFinished) && obj.iteration > 0;
        end
        
        function itIs = isTimeToRefresh(obj)
            itIs = mod(obj.iteration,obj.refreshInterval) == 0;
        end
        
        function createNamingManager(obj,cParams)
           nmS.costFuncNames = cParams.costFuncNames;
           nmS.costWeights = cParams.costWeights;
           nmS.constraintFuncs = cParams.constraintFuncs;
           nmS.optimizerName = cParams.optimizerName;
           
           obj.namingManager    = NamingManager(nmS);
        end
        
    end
    
    methods (Access = private, Static)
        
        function type = getChartType(title)
            if contains(title,'kappa')
                type = 'bar';
            elseif contains(title,'outit')
                type = 'stacked';
            elseif contains(title,'L2') || contains(title,'kkt')
                type = 'log';
            else
                type = 'plot';
            end
        end
        
    end
    
    methods (Access = ?Display_Abstract)
        
        function increaseRefreshIterator(obj)
            obj.iRefresh = obj.iRefresh + 1;
        end
        
    end
    
    methods
        
        function n = get.nFigs(obj)
            n = length(obj.figures);
        end
        
    end
    
end