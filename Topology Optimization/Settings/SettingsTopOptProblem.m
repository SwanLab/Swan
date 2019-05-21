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
        end
        
    end
    
    methods (Access = private)
        
        function setupProblemData(obj)
            s = obj.cParams.problemData;
            obj.problemData.problemFileName = s.problemFileName;
            obj.problemData.caseFileName = obj.loadedFile;
            obj.problemData.scale = s.scale;
            
            obj.createMesh();
            
            obj.problemData.pdim  = obj.mesh.pdim;
            obj.problemData.nelem = size(obj.mesh.connec,1);
        end
        
        function createMesh(obj)
            obj.mesh = Mesh_GiD(obj.problemData.problemFileName);
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
            obj.problemData.costFunctions       = obj.costSettings.getShapeFuncList();
            obj.problemData.costWeights         = obj.costSettings.weights();
            obj.problemData.constraintFunctions = obj.constraintSettings.getShapeFuncList();
            obj.problemData.nConstraints        = obj.constraintSettings.nShapeFuncs;
        end
        
        function createOptimizerSettings(obj)
            s = obj.cParams.optimizerSettings;
            obj.optimizerSettings = SettingsOptimizer(s);
            obj.optimizerSettings.problemData = obj.problemData;
            obj.optimizerSettings.initSettingsMonitorDocker(s);
            
            obj.optimizerSettings.nConstr = obj.problemData.nConstraints;
            obj.optimizerSettings.initSettingsHistoryPrinter(obj.problemData.caseFileName);
            obj.optimizerSettings.initSettingsPostProcess();
            
            obj.initOptimizerUnconstrainedSettings(s);
        end
        
        function initOptimizerUnconstrainedSettings(obj,cParams)
            s = cParams.uncOptimizerSettings;
            obj.optimizerSettings.uncOptimizerSettings = SettingsOptimizerUnconstrained(s);
            obj.optimizerSettings.uncOptimizerSettings.type = s.type;
            
            obj.initScalarProductSettings();
            obj.initLineSearchSettings(s);
        end
        
        function createVideoManagerSettings(obj)
            s.caseFileName  = obj.problemData.caseFileName;
            s.shallPrint    = obj.optimizerSettings.shallPrint;
            s.designVarType = obj.designVarSettings.type;
            s.pdim          = obj.problemData.pdim;
            obj.videoManagerSettings = SettingsVideoManager(s);
        end
        
        function initScalarProductSettings(obj)
            obj.optimizerSettings.uncOptimizerSettings.scalarProductSettings.filename = obj.problemData.problemFileName;
        end
        
        function initLineSearchSettings(obj,cParams)
            s = cParams.lineSearchSettings;
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings = SettingsLineSearch(s);
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.scalarProductSettings = obj.optimizerSettings.uncOptimizerSettings.scalarProductSettings;
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.optimizerType  = obj.optimizerSettings.uncOptimizerSettings.type;
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.filename       = obj.problemData.problemFileName;
        end
        
    end
    
end