classdef SettingsTopOptProblem < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsTopOptProblem.json'
    end
    
    properties (Access = public)
        settings
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
        isOld
    end
    
    methods (Access = public)
        
        function obj = SettingsTopOptProblem(varargin)
            if nargin == 1
                if ischar(varargin{1})
                    obj.loadParams(varargin{1});
                    cParams = obj.customParams;
                    obj.isOld = false;
                else
                    obj.settings = varargin{1};
                    cParams = [];
                    obj.isOld = true;
                end
            elseif nargin == 2
                obj.loadParams(varargin{1});
                cParams = obj.customParams;
                obj.settings = varargin{2};
                obj.isOld = obj.settings.isOld;
            end
            
            obj.setupProblemData(cParams);
            
            obj.createDesignVarSettings(cParams);
            obj.createHomogenizedVarComputerSettings();
            obj.createIncrementalSchemeSettings(cParams);
            obj.createCostSettings(cParams);
            obj.createConstraintSettings(cParams);
            obj.updateProblemData();
            obj.createOptimizerSettings(cParams);
            obj.createVideoManagerSettings();
        end
        
    end
    
    methods (Access = private)
        
        function setupProblemData(obj,cParams)
            if obj.isOld
                obj.problemData.problemFileName = obj.settings.filename;
                obj.problemData.caseFileName = obj.settings.case_file;
                obj.problemData.scale = obj.settings.ptype;
            else
                s = cParams.problemData;
                obj.problemData.problemFileName = s.problemFileName;
                obj.problemData.caseFileName = obj.loadedFile;
                obj.problemData.scale = s.scale;
            end
            
            obj.createMesh();
            
            obj.problemData.pdim  = obj.mesh.pdim;
            obj.problemData.nelem = size(obj.mesh.connec,1);
            if obj.isOld
                obj.settings.pdim = obj.problemData.pdim;
            end
        end
        
        function createMesh(obj)
            obj.mesh = Mesh_GiD(obj.problemData.problemFileName);
        end
        
        function createDesignVarSettings(obj,cParams)
            if obj.isOld
                s.type = obj.settings.designVariable;
                s.initialCase = obj.settings.initial_case;
                s.mesh = obj.mesh;
            else
                s = cParams.designVarSettings;
                s.mesh = obj.mesh;
            end
            obj.designVarSettings = SettingsDesignVariable(s);
            
            if obj.isOld
                obj.addGeomParams();
            end
        end
        
        function createHomogenizedVarComputerSettings(obj)
            if obj.isOld
                s.type  = obj.settings.homegenizedVariablesComputer;
                switch s.type
                    case 'ByInterpolation'
                        s.interpolation          = obj.settings.materialInterpolation;
                        s.typeOfMaterial         = obj.settings.material;
                        s.constitutiveProperties = obj.settings.TOL;
                    case 'ByVademecum'
                        s.fileName               = obj.settings.vademecumFileName;
                end
            else
                s = obj.homogenizedVarComputerSettings;
            end
            s.nelem = obj.problemData.nelem;
            s.dim   = obj.problemData.pdim;
            obj.homogenizedVarComputerSettings = SettingsHomogenizedVarComputer.create(s);
        end
        
        function createIncrementalSchemeSettings(obj,cParams)
            if obj.isOld
                s.nSteps = obj.settings.nsteps;
                s.shallPrintIncremental = obj.settings.printIncrementalIter;
                
                s.targetParamsSettings.VfracInitial = obj.settings.Vfrac_initial;
                s.targetParamsSettings.VfracFinal = obj.settings.Vfrac_final;
                s.targetParamsSettings.constrInitial = obj.settings.constr_initial;
                s.targetParamsSettings.constrFinal = obj.settings.constr_final;
                s.targetParamsSettings.optimalityInitial = obj.settings.optimality_initial;
                s.targetParamsSettings.optimalityFinal = obj.settings.optimality_final;
                s.targetParamsSettings.epsilonInitial = obj.settings.epsilon_initial;
                s.targetParamsSettings.epsilonFinal = obj.settings.epsilon_final;
                s.targetParamsSettings.epsilonIsoInitial = obj.settings.epsilon_isotropy_initial;
                s.targetParamsSettings.epsilonIsoFinal = obj.settings.epsilon_isotropy_final;
            else
                s = cParams.incrementalSchemeSettings;
            end
            s.mesh = obj.mesh;
            obj.incrementalSchemeSettings = SettingsIncrementalScheme(s);
        end
        
        function createCostSettings(obj,cParams)
            if obj.isOld
                s.weights = obj.settings.weights;
                for i = 1:length(obj.settings.cost)
                    sfS{i} = struct('type',obj.settings.cost{i});
                end
                s.shapeFuncSettings = sfS;
                s.problemData = obj.problemData;
                s.settings = obj.settings;
                obj.costSettings = SettingsCost_OLD(s);
            else
                s = cParams.costSettings;
                s.problemData = obj.problemData;
                obj.costSettings = SettingsCost(s);
            end
        end
        
        function createConstraintSettings(obj,cParams)
            if obj.isOld
                for i = 1:length(obj.settings.constraint)
                    sfS{i} = struct('type',obj.settings.constraint{i});
                end
                s.shapeFuncSettings = sfS;
                s.problemData = obj.problemData;
                s.settings = obj.settings;
                obj.constraintSettings = SettingsConstraint_OLD(s);
            else
                s = cParams.constraintSettings;
                s.problemData = obj.problemData;
                obj.constraintSettings = SettingsConstraint(s);
            end
        end
        
        function updateProblemData(obj)
            obj.problemData.costFunctions       = obj.costSettings.getShapeFuncList();
            obj.problemData.costWeights         = obj.costSettings.weights();
            obj.problemData.constraintFunctions = obj.constraintSettings.getShapeFuncList();
            obj.problemData.nConstraints        = obj.constraintSettings.nShapeFuncs;
        end
        
        function createOptimizerSettings(obj,cParams)
            obj.optimizerSettings = SettingsOptimizer();
            
            if obj.isOld
                s = [];
                obj.optimizerSettings.type = obj.settings.optimizer;
                obj.optimizerSettings.constraintCase = obj.settings.constraint_case;
                obj.optimizerSettings.maxIter = obj.settings.maxiter;
                obj.optimizerSettings.shallPrint = obj.settings.printing;
                obj.optimizerSettings.printMode  = obj.settings.printMode;
                obj.setupSettingsMonitor();
            else
                s = cParams.optimizerSettings;
                obj.optimizerSettings.problemData = obj.problemData;
                obj.optimizerSettings.type = s.type;
                obj.optimizerSettings.maxIter = s.maxIter;
                obj.optimizerSettings.shallPrint = s.shallPrint;
                obj.optimizerSettings.printMode  = s.printMode;
                obj.optimizerSettings.initSettingsMonitorDocker(s);
            end
            
            obj.optimizerSettings.nConstr = obj.problemData.nConstraints;
            obj.optimizerSettings.initSettingsHistoryPrinter(obj.problemData.caseFileName);
            obj.optimizerSettings.initSettingsPostProcess();
            
            obj.initOptimizerUnconstrainedSettings(s);
        end
        
        function initOptimizerUnconstrainedSettings(obj,cParams)
            if obj.isOld
                s = [];
                obj.optimizerSettings.uncOptimizerSettings.type = obj.settings.optimizerUnconstrained;
                obj.optimizerSettings.uncOptimizerSettings.ub = obj.settings.ub;
                obj.optimizerSettings.uncOptimizerSettings.lb = obj.settings.lb;
                obj.optimizerSettings.uncOptimizerSettings.e2 = obj.settings.e2;
            else
                s = cParams.uncOptimizerSettings;
                obj.optimizerSettings.uncOptimizerSettings = SettingsOptimizerUnconstrained(s);
                obj.optimizerSettings.uncOptimizerSettings.type = s.type;
            end
            
            obj.initScalarProductSettings();
            obj.initLineSearchSettings(s);
            
            obj.optimizerSettings.uncOptimizerSettings.filename            = obj.problemData.problemFileName;
            obj.optimizerSettings.uncOptimizerSettings.ptype               = obj.problemData.scale;
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
        
        function initLineSearchSettings(obj,cParams,spS)
            if obj.isOld
                obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.type            = obj.settings.line_search;
                obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.HJiter0         = obj.settings.HJiter0;
                obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.kfrac           = obj.settings.kfrac;
                obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.kappaMultiplier = obj.settings.kappaMultiplier;
            else
                s = cParams.lineSearchSettings;
                obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings = SettingsLineSearch(s);
            end
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.scalarProductSettings = obj.optimizerSettings.uncOptimizerSettings.scalarProductSettings;
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.optimizerType  = obj.optimizerSettings.uncOptimizerSettings.type;
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings.filename       = obj.problemData.problemFileName;
        end
        
        function setupSettingsMonitor(obj)
            obj.optimizerSettings.monitoringDockerSettings.showOptParams               = obj.settings.monitoring;
            obj.optimizerSettings.monitoringDockerSettings.refreshInterval             = obj.settings.monitoring_interval;
            obj.optimizerSettings.monitoringDockerSettings.shallDisplayDesignVar       = obj.settings.plotting;
            obj.optimizerSettings.monitoringDockerSettings.shallShowBoundaryConditions = obj.settings.showBC;
            
            obj.optimizerSettings.monitoringDockerSettings.optimizerName   = obj.settings.optimizer;
            obj.optimizerSettings.monitoringDockerSettings.problemID       = obj.settings.case_file;
            obj.optimizerSettings.monitoringDockerSettings.dim             = obj.settings.pdim;
            
            obj.optimizerSettings.monitoringDockerSettings.costFuncNames   = obj.problemData.costFunctions;
            obj.optimizerSettings.monitoringDockerSettings.costWeights     = obj.problemData.costWeights;
            obj.optimizerSettings.monitoringDockerSettings.constraintFuncs = obj.problemData.constraintFunctions;
        end
        
        function addGeomParams(obj)
            if ~isempty(obj.settings.levelSetDataBase)
                fields = fieldnames(obj.settings.levelSetDataBase);
                for i = 1:length(fields)
                    field = fields{i};
                    if isfield(obj.settings.levelSetDataBase,field)
                        obj.designVarSettings.levelSetCreatorSettings.(field) = obj.settings.levelSetDataBase.(field);
                    end
                end
            end
        end
        
    end
    
end