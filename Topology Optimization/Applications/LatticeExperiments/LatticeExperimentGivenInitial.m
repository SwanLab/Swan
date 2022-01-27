classdef LatticeExperimentGivenInitial < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
      m1
      m2
      alpha
    end
    
    properties (Access = private)
        fileNames
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = LatticeExperimentGivenInitial()
            obj.init();
            for icases = 1:numel(obj.fileNames)
                obj.readInitial();
                obj.createSettings(icases);
                obj.solveProblem();
                close all
            end
            
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileNames = {
                %'CantileverSymmetricFixingDirichletZone'
                %'CantileverSymmetricFixingMaxStressZone'
                'CantileverSymmetricWithoutFixing'};
        end
        
        function createSettings(obj,icases)
            s = SettingsTopOptProblem(obj.fileNames{icases});
            s.designVarSettings.creatorSettings.m1 = obj.m1;
            s.designVarSettings.creatorSettings.m2 = obj.m2;
            s.designVarSettings.creatorSettings.alpha0 = obj.alpha;
            obj.topOptSet = s;
        end
        
        function readInitial(obj)
            iteration = 262;
            fCase = 'CantileverSymmetricFixingMaxStressZone';
            folder = '/media/alex/My Passport/LatticeResults/CantileverNoEnglishFlag/CantileverSymmetricFixingMaxStressZone';
            s.fileName = [fCase,num2str(iteration)];
            s.folderPath = fullfile(folder);
            w = WrapperMshResFiles(s);
            w.compute();
            obj.m1 = w.dataRes.DesignVar1;
            obj.m2 = w.dataRes.DesignVar2;
            obj.alpha = w.dataRes.AlphaGauss';
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end
        
    end
    
end