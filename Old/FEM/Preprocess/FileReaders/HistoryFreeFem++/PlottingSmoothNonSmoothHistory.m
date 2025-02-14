classdef PlottingSmoothNonSmoothHistory < handle
    
    properties (Access = private)
        caseNames
        caseName
        iter
        color
        cost
        lambda
        resultsDir
    end
    
    methods (Access = public)
        
        function obj = PlottingSmoothNonSmoothHistory()
            obj.init();
            obj.readFiles();
            obj.plotHistory();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.resultsDir = '/home/alex/git-repos/OptimizationFreeFem/';
            obj.color = {'r','g'};
%            obj.caseNames = {'SmoothRectangle','Rectangle','RectangleWithMmax','VademecumSmoothCorner'};
            obj.caseNames = {'Rectangle','VademecumSmoothCorner'};
            
        end
        
        function obj = readFiles(obj)
            for icell = 1:numel(obj.caseNames)
                obj.iter = icell;
                obj.caseName = obj.caseNames{icell};
                obj.obtainCostAndLambda();
            end
        end
        
        function obtainCostAndLambda(obj)
            fullPath = [obj.resultsDir,[obj.caseName,'Results'],'/ctl/History.txt'];
            d.filePath = fullPath;
            rH = ReadingHistoryFile(d);
            obj.cost{obj.iter} = rH.cost;
            obj.lambda{obj.iter} = rH.lambda;
        end
        
        function plotHistory(obj)
            obj.plotVariable(obj.cost,'Cost');
            obj.plotVariable(obj.lambda,'Lambda');
        end
        
        function plotVariable(obj,var,name)
            f = figure;
            hold on
            for icell = 1:numel(var)
                h{icell} = plot(var{icell},[obj.color{icell},'+-']);
            end
            legend(obj.caseNames)
            pp = plotPrinter(f,h); 
            pp.print([obj.resultsDir,name])
        end
        
    end
    
end