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
                settings = varargin{1};
                cParams = [];
            elseif nargin == 2
                obj.loadParams(varargin{1});
                cParams = obj.customParams;
                settings = varargin{2};
            end
            obj.isOld = settings.isOld;
            obj.settings = settings;
            
            obj.setupProblemData(cParams);
            
            obj.createDesignVarSettings(cParams);
            obj.createHomogenizedVarComputerSettings();
            obj.createIncrementalSchemeSettings(cParams);
            obj.createCostSettings(cParams);
            obj.createConstraintSettings(cParams);
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
            obj.settings.pdim = obj.problemData.pdim;
            obj.problemData.costFunctions       = obj.settings.cost;
            obj.problemData.costWeights         = obj.settings.weights;
            obj.problemData.constraintFunctions = obj.settings.constraint;
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
            else
                s = cParams.incrementalSchemeSettings;
            end
            s.mesh = obj.mesh;
            s.settings = obj.settings;
            obj.incrementalSchemeSettings = SettingsIncrementalScheme(s);
        end
        
        function createCostSettings(obj,cParams)
            if obj.isOld
                s.weights = obj.settings.weights;
                for i = 1:length(obj.settings.cost)
                    sfS{i} = struct('type',obj.settings.cost{i});
                end
                s.shapeFuncSettings = sfS;
            else
                s = cParams.costSettings;
            end
            s.settings = obj.settings;
            s.problemData = obj.problemData;
            
            obj.costSettings = SettingsCost(s);
        end
        
        function createConstraintSettings(obj,cParams)
            if obj.isOld
                for i = 1:length(obj.settings.constraint)
                    sfS{i} = struct('type',obj.settings.constraint{i});
                end
                s.shapeFuncSettings = sfS;
            else
                s = cParams.constraintSettings;
            end
            s.settings = obj.settings;
            s.problemData = obj.problemData;
            obj.constraintSettings = SettingsConstraint(s);
        end
        
        function createOptimizerSettings(obj,cParams)
            obj.optimizerSettings = SettingsOptimizer();
            if obj.isOld
                s = [];
                obj.optimizerSettings.type = obj.settings.optimizer;
            else
                s = cParams.optimizerSettings;
                obj.optimizerSettings.type = s.type;
            end
            
            obj.optimizerSettings.nconstr              = obj.settings.nconstr;
            obj.optimizerSettings.target_parameters    = obj.settings.target_parameters;
            obj.optimizerSettings.constraint_case      = obj.settings.constraint_case;
            obj.optimizerSettings.maxiter              = obj.settings.maxiter;
            
            if obj.isOld
                obj.optimizerSettings.name             = obj.settings.optimizer;
                obj.optimizerSettings.shallPrint       = obj.settings.printing;
            else
                obj.optimizerSettings.name             = s.type;
                obj.optimizerSettings.shallPrint       = s.shallPrint;
            end
            obj.optimizerSettings.printMode            = obj.settings.printMode;
            
            obj.setupSettingsMonitor(s);
            obj.optimizerSettings.setupSettingsHistoryPrinter(obj.problemData.caseFileName);
            obj.optimizerSettings.setupSettingsPostProcess();
            
            obj.createOptimizerUnconstrainedSettings(s);
        end
        
        function createOptimizerUnconstrainedSettings(obj,cParams)
            obj.optimizerSettings.uncOptimizerSettings = SettingsOptimizerUnconstrained();
            if obj.isOld
                obj.optimizerSettings.uncOptimizerSettings.type = obj.settings.optimizerUnconstrained;
            else
                s = cParams.uncOptimizerSettings;
                obj.optimizerSettings.uncOptimizerSettings.type = s.type;
            end
            
            spS = obj.createScalarProductSettings();
            lsS = obj.createLineSearchSettings(spS);
            
            obj.optimizerSettings.uncOptimizerSettings.lineSearchSettings    = lsS;
            obj.optimizerSettings.uncOptimizerSettings.scalarProductSettings = spS;
            
            obj.optimizerSettings.uncOptimizerSettings.e2                  = obj.settings.e2;
            obj.optimizerSettings.uncOptimizerSettings.filter              = obj.settings.filter;
            obj.optimizerSettings.uncOptimizerSettings.printChangingFilter = obj.settings.printChangingFilter;
            obj.optimizerSettings.uncOptimizerSettings.filename            = obj.problemData.problemFileName;
            obj.optimizerSettings.uncOptimizerSettings.ptype               = obj.settings.ptype;
            obj.optimizerSettings.uncOptimizerSettings.lb                  = obj.settings.lb;
            obj.optimizerSettings.uncOptimizerSettings.ub                  = obj.settings.ub;
        end
        
        function createVideoManagerSettings(obj)
            s.caseFileName  = obj.problemData.caseFileName;
            s.shallPrint    = obj.optimizerSettings.shallPrint;
            s.designVarType = obj.designVarSettings.type;
            s.pdim          = obj.problemData.pdim;
            obj.videoManagerSettings = SettingsVideoManager(s);
        end
        
        function spS = createScalarProductSettings(obj)
            spS.filename = obj.problemData.problemFileName;
        end
        
        function lsS = createLineSearchSettings(obj,spS)
            lsS.scalarProductSettings = spS;
            lsS.line_search             = obj.settings.line_search;
            lsS.optimizerUnconstrained  = obj.optimizerSettings.uncOptimizerSettings.type;
            lsS.HJiter0         = obj.settings.HJiter0;
            lsS.filename        = obj.problemData.problemFileName;
            lsS.kappaMultiplier = obj.settings.kappaMultiplier;
            lsS.kfrac           = obj.settings.kfrac;
        end
        
        function setupSettingsMonitor(obj,cParams)
            if obj.isOld
                obj.optimizerSettings.settingsMonitor.showOptParams               = obj.settings.monitoring;
                obj.optimizerSettings.settingsMonitor.refreshInterval             = obj.settings.monitoring_interval;
                obj.optimizerSettings.settingsMonitor.shallDisplayDesignVar       = obj.settings.plotting;
                obj.optimizerSettings.settingsMonitor.shallShowBoundaryConditions = obj.settings.showBC;
                
                obj.optimizerSettings.settingsMonitor.optimizerName   = obj.settings.optimizer;
                obj.optimizerSettings.settingsMonitor.problemID       = obj.settings.case_file;
                obj.optimizerSettings.settingsMonitor.dim             = obj.settings.pdim;
            else
                s = cParams.settingsMonitor;
                obj.optimizerSettings.settingsMonitor = SettingsMonitoringDocker(s);
                
                obj.optimizerSettings.settingsMonitor.optimizerName   = obj.optimizerSettings.name;
                obj.optimizerSettings.settingsMonitor.problemID       = obj.problemData.caseFileName;
                obj.optimizerSettings.settingsMonitor.dim             = obj.problemData.pdim;
            end
            
            obj.optimizerSettings.settingsMonitor.costFuncNames   = obj.problemData.costFunctions;
            obj.optimizerSettings.settingsMonitor.costWeights     = obj.problemData.costWeights;
            obj.optimizerSettings.settingsMonitor.constraintFuncs = obj.problemData.constraintFunctions;
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