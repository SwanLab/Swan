classdef PlottingSmoothNonSmoothHistory < handle
    
    properties (Access = private)
        costSmooth
        lambdaSmooth
        costNonSmooth
        lambdaNonSmooth
        costNonSmoothMmax
        lambdaNonSmoothMmax
        resultsDir        
    end
    
    methods (Access = public)
        
        function obj = PlottingSmoothNonSmoothHistory()
            obj.init();
            obj.readSmooth();
            obj.readNonSmooth();
            obj.readNonSmoothMmax();            
            obj.plotHistory();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.resultsDir = '/home/alex/Desktop/OptimizationFreeFem/';
        end
        
        function obj = readSmooth(obj)
            caseName = 'SmoothRectangleResults';
            fullPath = [obj.resultsDir,caseName,'/ctl/History.txt'];
            d.filePath = fullPath;
            rH = ReadingHistoryFile(d);
            obj.costSmooth = rH.cost;
            obj.lambdaSmooth = rH.lambda;
        end
        
        function obj = readNonSmooth(obj)
            caseName = 'RectangleResults';
            fullPath = [obj.resultsDir,caseName,'/ctl/History.txt'];
            d.filePath = fullPath;
            rH = ReadingHistoryFile(d);
            obj.costNonSmooth = rH.cost;
            obj.lambdaNonSmooth = rH.lambda;
        end
        
        function obj = readNonSmoothMmax(obj)
            caseName = 'RectangleWithMmaxResults';
            fullPath = [obj.resultsDir,caseName,'/ctl/History.txt'];
            d.filePath = fullPath;
            rH = ReadingHistoryFile(d);
            obj.costNonSmoothMmax = rH.cost;
            obj.lambdaNonSmoothMmax = rH.lambda;
        end
        
        
        
        
        function plotHistory(obj)
            f = figure;
            hold on
            h{1} = plot(obj.costSmooth,'b+-');
            h{2} = plot(obj.costNonSmooth,'r+-');
            h{3} = plot(obj.costNonSmoothMmax,'g+-');
            legend('SmoothCost','NonSmoothCost','NonSmoothCostMmax')
            
            pp = plotPrinter(f,h); 
            pp.print([obj.resultsDir,'Cost'])
            
            f = figure;
            hold on            
            h{1} = plot(obj.lambdaSmooth,'b+-');
            h{2} = plot(obj.lambdaNonSmooth,'r+-');
            h{3} = plot(obj.lambdaNonSmoothMmax,'g+-');            
            legend('SmoothLambda','NonSmoothLambda','NonSmoothCostMmax')   
            
            pp = plotPrinter(f,h); 
            pp.print([obj.resultsDir,'Lambda'])            
        end
        
    end
    
    
    
end