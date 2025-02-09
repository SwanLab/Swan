classdef PerimeterTopOptExperiment < handle

    properties (Access = private)
        fileName
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = PerimeterTopOptExperiment()
            obj.init();              
            obj.createTopOptSettings();            
            obj.solveProblem();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'InteriorPerimeter';
        end
        
        function createTopOptSettings(obj)
          t = SettingsTopOptProblem(obj.fileName);            
          ls = obj.readLevelSet();
         
          t.designVarSettings.initialCase = 'given';
          t.designVarSettings.creatorSettings.value = ls;               
          vF = t.incrementalSchemeSettings.targetParamsSettings.VfracFinal;            
          t.incrementalSchemeSettings.targetParamsSettings.VfracInitial = vF;            
          obj.topOptSet = t;              
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end        
        
        function ls = readLevelSet(obj)
            iter = 200;
            fCase = 'InteriorPerimeter';
            folder = '/media/alex/My Passport/PerimeterResults/FineMesh/InitialTopology/';            
            s.fileName = [fCase,num2str(iter)];
            s.folderPath = fullfile(folder);            
            wM = WrapperMshResFiles(s);
            wM.compute();
            ls = wM.dataRes.DesignVar1;
        end
        
    end
    
end