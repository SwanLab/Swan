classdef SettingsTopOptProblem < AbstractSettings
    
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
            elseif nargin == 2
                obj.loadParams(varargin{1});
                settings = varargin{2};
            end
            obj.isOld = settings.isOld;
            obj.settings = settings;
            
            obj.setupProblemData();
            
            obj.createDesignVarSettings();
            obj.createHomogenizedVarComputerSettings();
            obj.createIncrementalSchemeSettings();
            obj.createCostSettings();
            obj.createConstraintSettings();
            obj.createOptimizerSettings();
            obj.createVideoManagerSettings();
        end
        
    end
    
    methods (Access = private)
        
        function createMesh(obj)
            obj.mesh = Mesh_GiD(obj.problemData.problemFileName);
        end
        
        function createDesignVarSettings(obj)
            obj.designVarSettings.mesh = obj.mesh;
            if obj.isOld
                obj.designVarSettings.type = obj.settings.designVariable;
                obj.designVarSettings.initialCase = obj.settings.initial_case;
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
        
        function setupProblemData(obj)
            if obj.isOld
                obj.problemData.problemFileName = obj.settings.filename;
                obj.problemData.caseFileName = obj.settings.case_file;
                obj.problemData.scale = obj.settings.ptype;
            else
                obj.problemData.caseFileName = obj.loadedFile;
            end
            
            obj.createMesh();
            
            obj.problemData.pdim  = obj.mesh.pdim;
            obj.problemData.nelem = size(obj.mesh.connec,1);
            obj.settings.pdim = obj.problemData.pdim;
            obj.problemData.costFunctions       = obj.settings.cost;
            obj.problemData.costWeights         = obj.settings.weights;
            obj.problemData.constraintFunctions = obj.settings.constraint;
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
        
        function createIncrementalSchemeSettings(obj)
            tpS = obj.createTargetParamsSettings();
            obj.incrementalSchemeSettings.settingsTargetParams = tpS;
            if obj.isOld
                obj.incrementalSchemeSettings.nSteps = obj.settings.nsteps;
            end
            obj.incrementalSchemeSettings.shallPrintIncremental = obj.settings.printIncrementalIter;
            
            obj.incrementalSchemeSettings.mesh = obj.mesh;
        end
        
        function createCostSettings(obj)
            if obj.isOld
                weights = obj.settings.weights;
                for i = 1:length(obj.settings.cost)
                    sfS{i} = struct('type',obj.settings.cost{i});
                end
            else
                sfS = obj.costSettings.shapeFuncSettings;
                weights = obj.costSettings.weights;
            end
            obj.costSettings = struct;
            obj.costSettings.settings = obj.settings;
            obj.costSettings.weights  = weights;
            obj.costSettings.shapeFuncSettings = obj.createShapeFunctionsSettings(sfS);
            obj.costSettings.nShapeFuncs = length(obj.costSettings.shapeFuncSettings);
            
        end
        
        function createConstraintSettings(obj)
            if obj.isOld
                for i = 1:length(obj.settings.constraint)
                    sfS{i} = struct('type',obj.settings.constraint{i});
                end
            else
                sfS = obj.constraintSettings.shapeFuncSettings;
            end
            obj.constraintSettings = struct;
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
        
        function createOptimizerSettings(obj)
            obj.settings.pdim = obj.problemData.pdim;
            
            uoS = obj.createOptimizerUnconstrainedSettings();
            obj.optimizerSettings.uncOptimizerSettings = uoS;
            
            if obj.isOld
                obj.optimizerSettings.type = obj.settings.optimizer;
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
            
            obj.setupSettingsMonitor();
            obj.optimizerSettings.setupSettingsHistoryPrinter(obj.problemData.caseFileName);
            obj.optimizerSettings.setupSettingsPostProcess();
        end
        
        function uoS = createOptimizerUnconstrainedSettings(obj)
            uoS = obj.optimizerSettings.uncOptimizerSettings;
            if obj.isOld
                uoS.type = obj.settings.optimizerUnconstrained;
            end
            
            spS = obj.createScalarProductSettings();
            lsS = obj.createLineSearchSettings(spS);
            
            uoS.lineSearchSettings    = lsS;
            uoS.scalarProductSettings = spS;
            
            uoS.e2                  = obj.settings.e2;
            uoS.filter              = obj.settings.filter;
            uoS.printChangingFilter = obj.settings.printChangingFilter;
            uoS.filename            = obj.problemData.problemFileName;
            uoS.ptype               = obj.settings.ptype;
            uoS.lb                  = obj.settings.lb;
            uoS.ub                  = obj.settings.ub;
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
        
        function setupSettingsMonitor(obj)
            if obj.isOld
                obj.optimizerSettings.settingsMonitor.showOptParams               = obj.settings.monitoring;
                obj.optimizerSettings.settingsMonitor.refreshInterval             = obj.settings.monitoring_interval;
                obj.optimizerSettings.settingsMonitor.shallDisplayDesignVar       = obj.settings.plotting;
                obj.optimizerSettings.settingsMonitor.shallShowBoundaryConditions = obj.settings.showBC;
                
                obj.optimizerSettings.settingsMonitor.optimizerName   = obj.settings.optimizer;
                obj.optimizerSettings.settingsMonitor.problemID       = obj.settings.case_file;
                obj.optimizerSettings.settingsMonitor.dim             = obj.settings.pdim;
            else
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