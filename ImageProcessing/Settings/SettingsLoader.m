classdef SettingsLoader < handle
    
    properties (Access = public)
        settings
        
    end
    
    properties (Access = private)
        possibleInput
    end
    
    
    methods (Access = public)
        
        function obj = SettingsLoader(problemName)
            obj.loadListOfPossibleInputData();
            obj.loadInputData(problemName);
        end
        
        function loadListOfPossibleInputData(obj)
            input = {'imageFile',...
                'lipschitzConstant',...
                'totalVariationWeigth',...
                'noiseAmplitud',...
                'maxIter',...
                'optimizer'};
            obj.possibleInput = input;
        end
        
        function loadInputData(obj,problemName)
            run(problemName);
            allData = who;
            for ivar = 1:length(allData)
                data = allData{ivar};
                s.(data) = eval(data);
            end
            obj.settings = s;
        end
        
    end
    
end