classdef ConvergencePlotComputerMaterialDesignExperimentCircle < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        filesPath
        testNames
        linesData
        barData
        fieldsData
        outPutPlotPath
        alphas
        legends
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConvergencePlotComputerMaterialDesignExperimentCircle()
            obj.init();
            obj.loadFieldsData();
            obj.plotData();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.filesPath = '/media/alex/MyPassport/MaterialDesign/CStar/';
            obj.outPutPlotPath = '/home/alex/Dropbox/MaterialDesign/CC';
            obj.linesData = [1:6,9:15];
            obj.barData = [7:8];
            obj.alphas = {'0.10','0.00'};
            obj.legends = {'$\alpha=0.1$','$\alpha=0$'};
            obj.testNames = {'LargeCircleWithPerimeter10_1b/','LargeCircleWithPerimeter10_12b'};
        end
        
        function loadFieldsData(obj)
            for iCase = 1:numel(obj.testNames)
                testPath = fullfile(obj.filesPath,obj.testNames{iCase},'/');
                s.testPath  = testPath;
                s.linesData = obj.linesData;
                s.barData   = obj.barData;
                m = MonitoringDataLoader(s);
                obj.fieldsData{iCase} = m.obtainData();
            end
        end
        
        function plotData(obj)
            obj.plotCost();
            obj.plotPerimeter();
            obj.plotVolum();
            obj.plotEpsilon();
        end
        
        function plotCost(obj)
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                fieldToPlot = 'C - C not scaled';
                [x,y] = obj.obtainField(fieldToPlot,fieldData);
                p{iCase} = semilogy(x,y);
                hold on
            end
            legend(obj.legends,'Interpreter','latex','Location','Best');
            f = figure(1);
            p = plotPrinter(f,p);
            p.print(fullfile(obj.outPutPlotPath,'Cost'))
        end
        
        function plotPerimeter(obj)
            figure(2)
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                fieldToPlot = 'Perimeter non scaled';
                %fieldToPlot = ['PerimeterInterior (wt. ',obj.alphas{iCase},')'];
                [x,y] = obj.obtainField(fieldToPlot,fieldData);
                p{iCase} = plot(x,y);
                hold on
            end               
            legend(obj.legends,'Interpreter','latex','Location','Best');
            f = figure(2);
            p = plotPrinter(f,p);
            p.print(fullfile(obj.outPutPlotPath,'Perimeter'))
        end      
        
        function plotVolum(obj)
            figure(3)
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                fieldToPlot = 'Volum';
                [x,y] = obj.obtainField(fieldToPlot,fieldData);
                p{iCase} = plot(x,y);
                hold on
            end   
            legend(obj.legends,'Interpreter','latex','Location','Best');
            f = figure(3);
            p = plotPrinter(f,p);
            p.print(fullfile(obj.outPutPlotPath,'Volume'))
        end            
        
        function plotEpsilon(obj)
            figure(4)
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                [x,y] = obj.obtainField('epsilon over h',fieldData);
                p{iCase} = plot(x,y);
                hold on
            end
            legend(obj.legends,'Interpreter','latex','Location','Best');
            f = figure(4);
            p = plotPrinter(f,p);
            p.print(fullfile(obj.outPutPlotPath,'Epsilon'))
        end
        
    end
    
    methods (Access = private, Static)

        function [xV,yV] = obtainField(fieldName,fieldData)
            for iField = 1:numel(fieldData)
                title = fieldData{iField}.title;
                if strcmp(fieldName,title)
                    xV = fieldData{iField}.xValue;
                    yV = fieldData{iField}.yValue;
                end
            end
        end
        
    end
    
end