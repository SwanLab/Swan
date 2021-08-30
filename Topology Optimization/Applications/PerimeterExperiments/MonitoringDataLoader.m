classdef MonitoringDataLoader < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        filesPath
        testName
    end
    
    properties (Access = private)        
        testPath
        subFiguresData
        subTitleData        
        figureTitle
        figureData
        linesData
        barData
    end
    
    methods (Access = public)
        
        function obj = MonitoringDataLoader(cParams)
            obj.init(cParams);
        end
        
        function fieldsD = obtainData(obj)
            obj.loadSubFiguresDataFromMonitoring();
            obj.loadSubFiguresTitlesFromMonitoring();
            nSubFigures = length(obj.subFiguresData);
            fieldsD = cell(nSubFigures,1);
            for isubF = 1:nSubFigures
                fieldsD{isubF}.xValue = obj.loadXvalue(isubF);
                fieldsD{isubF}.yValue = obj.loadYvalue(isubF);
                fieldsD{isubF}.title  = obj.loadTitle(isubF);
                fieldsD{isubF}.title
            end
            close all;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.testPath  = cParams.testPath;
            obj.linesData = cParams.linesData;
            obj.barData   = cParams.barData;
        end
        
        function loadSubFiguresDataFromMonitoring(obj)
            obj.loadLinesData();     
            obj.loadBarData();
        end
        
        function loadLinesData(obj)
            obj.loadData(obj.linesData,'line')
        end
        
        function loadBarData(obj)
            obj.loadData(obj.barData,'bar')
        end        
        
        function loadData(obj,lines,name)
            fNameMon = fullfile(obj.testPath,'Monitoring.fig');
            obj.figureData = openfig(fNameMon,'invisible');                        
            lineData  = findobj(obj.figureData,'Type',name);              
            for i = 1:numel(lines)
                obj.subFiguresData{lines(i)} = lineData(i);
            end 
        end
        
        function loadSubFiguresTitlesFromMonitoring(obj)
            fNameMon = fullfile(obj.testPath,'Monitoring.fig');
            obj.figureTitle = openfig(fNameMon,'invisible');           
            subFig = get(obj.figureTitle,'children');
            obj.subTitleData = get(subFig,'Title');
        end
        
        function field = loadXvalue(obj,number)
            subF = obj.subFiguresData{number};
            field = get(subF,'Xdata');            
        end
        
        function field = loadYvalue(obj,number)
            subF = obj.subFiguresData{number};            
            field = get(subF,'Ydata');            
        end        
        
        function title = loadTitle(obj,number)
            subF = obj.subTitleData{number};            
            title = subF.String;            
        end                
        
    end
    
end