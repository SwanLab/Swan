classdef SettingsTopOptProblem < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTopOptProblem.json'
    end
    
    properties (Access = public)
        fileName
        problemData = struct
        designVarSettings
        homogenizedVarComputerSettings
        incrementalSchemeSettings
        costSettings
        constraintSettings
        optimizerSettings
        videoMakerSettings
    end
    
    methods (Access = public)
        
        function obj = SettingsTopOptProblem(fileName)
            settings = Settings(fileName);
            translator = SettingsTranslator();
            translator.translate(settings);
            fileName = translator.fileName;
            if nargin == 1
                obj.loadParams([fileName,'.json']);
                obj.fileName = fileName;
            end
            obj.setupProblemData();
            obj.createDesignVarSettings();
            obj.createHomogenizedVarComputerSettings();
            obj.createIncrementalSchemeSettings();
            obj.createCostSettings();
            obj.createConstraintSettings();
            obj.updateProblemData();
            obj.createOptimizerSettings();
            obj.createVideoManagerSettings();
            %obj.printSummary();
        end
        
    end
    
    methods (Access = private)
        
        function setupProblemData(obj)
            s = obj.cParams.problemData;
            s.caseFileName = obj.fileName;
            obj.problemData = TopOptProblemDataContainer(s);
        end
        
        function createDesignVarSettings(obj)
            s = obj.cParams.designVarSettings;
            s.femData = obj.problemData.femData;
            obj.designVarSettings = SettingsDesignVariable(s);
        end
        
        function createHomogenizedVarComputerSettings(obj)
            s = obj.homogenizedVarComputerSettings;
            s.nelem = obj.problemData.femData.nelem;
            s.dim   = obj.problemData.femData.dim;
            obj.homogenizedVarComputerSettings = SettingsHomogenizedVarComputer.create(s);
        end
        
        function createIncrementalSchemeSettings(obj)
            s = obj.cParams.incrementalSchemeSettings;
            obj.incrementalSchemeSettings = SettingsIncrementalScheme(s);
        end
        
        function createCostSettings(obj)
            s = obj.cParams.costSettings;
            s.femData = obj.problemData.femData;
            obj.costSettings = SettingsCost(s);
        end
        
        function createConstraintSettings(obj)
            s = obj.cParams.constraintSettings;
            s.femData = obj.problemData.femData;
            obj.constraintSettings = SettingsConstraint(s);
        end
        
        function updateProblemData(obj)
            s.costFunctions       = obj.costSettings.getShapeFuncList();
            s.costWeights         = obj.costSettings.weights();
            s.constraintFunctions = obj.constraintSettings.getShapeFuncList();
            s.nConstraints        = obj.constraintSettings.nShapeFuncs;
            obj.problemData.loadParams(s);
        end
        
        function createOptimizerSettings(obj)
            s = obj.cParams.optimizerSettings;
            obj.optimizerSettings = SettingsOptimizer(s);
            
            s2.problemData = obj.problemData;
            s2.nConstr = obj.problemData.nConstraints;
            obj.optimizerSettings.loadParams(s2);
            obj.optimizerSettings.init();
        end
        
        function createVideoManagerSettings(obj)
            s.shallPrint    = obj.optimizerSettings.shallPrint;
            s.designVarType = obj.designVarSettings.type;
            s.pdim          = obj.problemData.femData.dim;
            s.caseFileName  = obj.fileName;            
            obj.videoMakerSettings = SettingsVideoMaker(s);
        end
        
        function printSummary(obj)
            if obj.isNotTest()
                fprintf('<strong>%s</strong>\n\n',obj.problemData.caseFileName)
                fprintf('\t-Optimizer: <strong>%s</strong>\n',obj.optimizerSettings.type);
                if strcmp(obj.optimizerSettings.type,'AlternatingPrimalDual')
                    fprintf('\t-Primal Updater: <strong>%s</strong>\n',obj.optimizerSettings.uncOptimizerSettings.type);
                end
                fprintf('\t-Cost: <strong>%s</strong>, ',obj.problemData.costFunctions{:})
                fprintf('\n\t-Constraints: ')
                fprintf('<strong>%s</strong>, ', obj.problemData.constraintFunctions{:})
                fprintf('\n\t-Incremental Steps: <strong>%.0f</strong> \n ',obj.incrementalSchemeSettings.nSteps)
                fprintf('\t-Max Iters: <strong>%.0f</strong> \n ',obj.optimizerSettings.maxIter)
            end
        end
        
        function itIsNot = isNotTest(obj)
            itIsNot = ~contains(obj.problemData.caseFileName,'test','IgnoreCase',true);
        end
        
    end
    
end