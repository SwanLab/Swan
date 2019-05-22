classdef SettingsTopOptProblem < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTopOptProblem.json'
    end
    
    properties (Access = public)
        problemData = struct
        designVarSettings
        homogenizedVarComputerSettings
        incrementalSchemeSettings
        costSettings
        constraintSettings
        optimizerSettings
        videoManagerSettings
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsTopOptProblem(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
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
            obj.printSummary();
        end
        
    end
    
    methods (Access = private)
        
        function setupProblemData(obj)
            s = obj.cParams.problemData;
            obj.problemData = TopOptProblemDataContainer(s);
            obj.problemData.caseFileName = obj.loadedFile;
            obj.problemData.femFileName = s.femFileName;
            obj.problemData.scale = s.scale;
            
            obj.createMesh();
            
            obj.problemData.pdim  = obj.mesh.pdim;
            obj.problemData.nelem = size(obj.mesh.connec,1);
        end
        
        function createMesh(obj)
            obj.mesh = Mesh_GiD(obj.problemData.femFileName);
        end
        
        function createDesignVarSettings(obj)
            s = obj.cParams.designVarSettings;
            s.mesh = obj.mesh;
            obj.designVarSettings = SettingsDesignVariable(s);
        end
        
        function createHomogenizedVarComputerSettings(obj)
            s = obj.homogenizedVarComputerSettings;
            s.nelem = obj.problemData.nelem;
            s.dim   = obj.problemData.pdim;
            obj.homogenizedVarComputerSettings = SettingsHomogenizedVarComputer.create(s);
        end
        
        function createIncrementalSchemeSettings(obj)
            s = obj.cParams.incrementalSchemeSettings;
            s.mesh = obj.mesh;
            obj.incrementalSchemeSettings = SettingsIncrementalScheme(s);
        end
        
        function createCostSettings(obj)
            s = obj.cParams.costSettings;
            s.problemData = obj.problemData;
            obj.costSettings = SettingsCost(s);
        end
        
        function createConstraintSettings(obj)
            s = obj.cParams.constraintSettings;
            s.problemData = obj.problemData;
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
            s.caseFileName  = obj.problemData.caseFileName;
            s.shallPrint    = obj.optimizerSettings.shallPrint;
            s.designVarType = obj.designVarSettings.type;
            s.pdim          = obj.problemData.pdim;
            obj.videoManagerSettings = SettingsVideoManager(s);
        end
        
        function printSummary(obj)
            if obj.isNotTest()
                fprintf('<strong>%s</strong>\n',obj.problemData.caseFileName)
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