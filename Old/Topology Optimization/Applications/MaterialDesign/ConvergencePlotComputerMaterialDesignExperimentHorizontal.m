classdef ConvergencePlotComputerMaterialDesignExperimentHorizontal < handle
    
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
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConvergencePlotComputerMaterialDesignExperimentHorizontal()
            obj.init();
            obj.loadFieldsData();
            obj.plotData();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.filesPath = '/media/alex/MyPassport/MaterialDesign/CStar/';
            obj.outPutPlotPath = '/home/alex/Dropbox/MaterialDesign/CC/Horizontal/';
            obj.linesData = [1:2,5:10];
            obj.barData = [3:4];
            obj.testNames = {'HorizontalFromCircleWithPerimeter10_8_LargeCircle/'};
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
            p{1} = obj.plotCost();
            p{2} = obj.plotVolum();
            f = figure(1);
            legend({'$||\bf{C}(\rho) - {C}^*||_2$','$\textrm{Vol}(\rho)$'},'Interpreter','latex','Location','Best');
            p = plotPrinter(f,p);
            p.print(fullfile(obj.outPutPlotPath,'CostAndVolume'))            
        end
        
        function p = plotCost(obj)
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                fieldToPlot = 'C - C not scaled';
                [x,y] = obj.obtainField(fieldToPlot,fieldData);
                p = semilogy(x,y);
                hold on
            end              
        end
        
        function p = plotVolum(obj)
            figure(1)
            hold on
            yyaxis right
            for iCase = 1:numel(obj.testNames)
                fieldData = obj.fieldsData{iCase};
                fieldToPlot = 'Volum';
                [x,y] = obj.obtainField(fieldToPlot,fieldData);
                p = plot(x,y);
                ylim([0,1])
                hold on
            end   
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