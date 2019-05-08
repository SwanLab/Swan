classdef SettingsTopOptProblem < AbstractSettings_B
    
    properties (Access = protected)
        defaultParamsName = 'paramsTopOptProblem'
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
                %                 obj.loadParams(varargin{1});
                cParams  = varargin{1};
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
                obj.problemData.scale = s.problemData.scale;
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
            obj.designVarSettings = SettingsDesignVariable();
            if obj.isOld
                obj.designVarSettings.type = obj.settings.designVariable;
                obj.designVarSettings.initialCase = obj.settings.initial_case;
                obj.designVarSettings.mesh = obj.mesh;
            else
                s = cParams.designVarSettings;
                s.mesh = obj.mesh;
                obj.designVarSettings.type = s.type;
                obj.designVarSettings.initialCase = s.initialCase;
                obj.designVarSettings.mesh = s.mesh;
            end
            obj.designVarSettings.levelSetCreatorSettings       = obj.settings.levelSetDataBase;
            obj.designVarSettings.levelSetCreatorSettings.ndim  = obj.mesh.ndim;
            obj.designVarSettings.levelSetCreatorSettings.coord = obj.mesh.coord;
            obj.designVarSettings.levelSetCreatorSettings.type = obj.designVarSettings.initialCase;
            switch obj.designVarSettings.levelSetCreatorSettings.type
                case 'holes'
                    obj.designVarSettings.levelSetCreatorSettings.dirichlet = obj.mesh.dirichlet;
                    obj.designVarSettings.levelSetCreatorSettings.pointload = obj.mesh.pointload;
            end
        end
        
        function createHomogenizedVarComputerSettings(obj)
            obj.homogenizedVarComputerSettings.type                   = obj.settings.homegenizedVariablesComputer;
            obj.homogenizedVarComputerSettings.interpolation          = obj.settings.materialInterpolation;
            obj.homogenizedVarComputerSettings.typeOfMaterial         = obj.settings.material;
            obj.homogenizedVarComputerSettings.constitutiveProperties = obj.settings.TOL;
            obj.homogenizedVarComputerSettings.vademecumFileName      = obj.settings.vademecumFileName;
            obj.homogenizedVarComputerSettings.dim                    = obj.problemData.pdim;
            obj.homogenizedVarComputerSettings.nelem                  = obj.problemData.nelem;
        end
        
        function createIncrementalSchemeSettings(obj,cParams)
            tpS = obj.createTargetParamsSettings();
            obj.incrementalSchemeSettings = SettingsIncrementalScheme();
            obj.incrementalSchemeSettings.settingsTargetParams = tpS;
            if obj.isOld
                obj.incrementalSchemeSettings.nSteps = obj.settings.nsteps;
            else
                s = cParams.incrementalSchemeSettings;
                obj.incrementalSchemeSettings.nSteps = s.nSteps;
            end
            obj.incrementalSchemeSettings.shallPrintIncremental = obj.settings.printIncrementalIter;
            obj.incrementalSchemeSettings.mesh = obj.mesh;
        end
        
        function createCostSettings(obj,cParams)
            obj.costSettings = SettingsCost();
            if obj.isOld
                obj.costSettings.weights = obj.settings.weights;
                for i = 1:length(obj.settings.cost)
                    sfS{i} = struct('type',obj.settings.cost{i});
                end
            else
                s = cParams.costSettings;
                obj.costSettings.weights = s.weights;
                sfS = obj.costSettings.shapeFuncSettings;
            end
            obj.costSettings.settings = obj.settings;
            obj.costSettings.shapeFuncSettings = obj.createShapeFunctionsSettings(sfS);
            obj.costSettings.nShapeFuncs = length(obj.costSettings.shapeFuncSettings);
        end
        
        function createConstraintSettings(obj,cParams)
            obj.constraintSettings = SettingsConstraint();
            if obj.isOld
                for i = 1:length(obj.settings.constraint)
                    sfS{i} = struct('type',obj.settings.constraint{i});
                end
            else
                s = cParams.constraintSettings;
                sfS = s.shapeFuncSettings;
            end
            obj.constraintSettings.settings = obj.settings;
            obj.constraintSettings.shapeFuncSettings = obj.createShapeFunctionsSettings(sfS);
            obj.constraintSettings.nShapeFuncs = length(obj.constraintSettings.shapeFuncSettings);
        end
        
        function cParams = createShapeFunctionsSettings(obj,s)
            nSF = length(s);
            cParams = cell(nSF,1);
            for iSF = 1:nSF
                s{iSF}.filename = obj.problemData.problemFileName;
                s{iSF}.scale = obj.problemData.scale;
                s{iSF}.filterParams = obj.createFilterSettings();
                cParams{iSF} = SettingsShapeFunctional().create(s{iSF},obj.settings);
            end
        end
        
        function s = createFilterSettings(obj)
            s = SettingsFilter();
            s.filterType = obj.settings.filter;
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
            obj.videoManagerSettings.caseFileName  = obj.problemData.caseFileName;
            obj.videoManagerSettings.shallPrint    = obj.optimizerSettings.shallPrint;
            obj.videoManagerSettings.designVarType = obj.designVarSettings.type;
            obj.videoManagerSettings.pdim          = obj.problemData.pdim;
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
                obj.optimizerSettings.settingsMonitor.showOptParams               = s.monitoring;
                obj.optimizerSettings.settingsMonitor.refreshInterval             = s.monitoring_interval;
                obj.optimizerSettings.settingsMonitor.shallDisplayDesignVar       = s.plotting;
                obj.optimizerSettings.settingsMonitor.shallShowBoundaryConditions = s.showBC;
                
                obj.optimizerSettings.settingsMonitor.optimizerName   = obj.optimizerSettings.name;
                obj.optimizerSettings.settingsMonitor.problemID       = obj.problemData.caseFileName;
                obj.optimizerSettings.settingsMonitor.dim             = obj.problemData.pdim;
            end
            
            obj.optimizerSettings.settingsMonitor.costFuncNames   = obj.problemData.costFunctions;
            obj.optimizerSettings.settingsMonitor.costWeights     = obj.problemData.costWeights;
            obj.optimizerSettings.settingsMonitor.constraintFuncs = obj.problemData.constraintFunctions;
        end
        
        function tpS = createTargetParamsSettings(obj)
            tpS = SettingsTargetParamsManager;
            tpS.VfracInitial = obj.settings.Vfrac_initial;
            tpS.VfracFinal = obj.settings.Vfrac_final;
            tpS.constrInitial = obj.settings.constr_initial;
            tpS.constrFinal = obj.settings.constr_final;
            tpS.optimalityInitial = obj.settings.optimality_initial;
            tpS.optimalityFinal = obj.settings.optimality_final;
            tpS.epsilonInitial = obj.settings.epsilon_initial;
            tpS.epsilonFinal = obj.settings.epsilon_final;
            tpS.epsilonIsotropyInitial = obj.settings.epsilon_isotropy_initial;
            tpS.epsilonIsotropyFinal = obj.settings.epsilon_isotropy_final;
        end
        
    end
    
end