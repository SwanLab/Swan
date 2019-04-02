classdef IncrementalScheme < handle
    
    properties (GetAccess = public, SetAccess = private)
        iStep
        nSteps
        targetParams
    end
    
    properties (Access = private)
        targetParamsManager
        
        epsilonInitial
        epsilonFinal
        epsilonPerInitial
        epsilonPerFinal
        
        shallDisplayStep
    end
    
    methods (Access = public)
        
        function obj = IncrementalScheme(settings,mesh)
            obj.init(settings,mesh);
            obj.createTargetParams(settings);
        end
        
        function next(obj)
            obj.incrementStep();
            obj.updateTargetParams();
        end
        
        
        function display(obj)
            disp(['Incremental Scheme - Step: ',int2str(obj.iStep),' of ',int2str(obj.nSteps)]);
        end
        
        function itDoes = hasNext(obj)
            if obj.iStep < obj.nSteps
                itDoes = true;
            else
                itDoes = false;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,settings,mesh)
            obj.iStep = 0;
            obj.nSteps = settings.nsteps;
            obj.setupEpsilons(settings,mesh);
            obj.setWhetherShallDisplayStep(settings);
        end
        
        function createTargetParams(obj,settings)
            settingsTargetParams = struct;
            settingsTargetParams.nSteps = obj.nSteps;
            settingsTargetParams.scale = settings.ptype;
            settingsTargetParams.Vfrac_initial = settings.Vfrac_initial;
            settingsTargetParams.Vfrac_final = settings.Vfrac_final;
            settingsTargetParams.constr_initial = settings.constr_initial;
            settingsTargetParams.constr_final = settings.constr_final;
            settingsTargetParams.optimality_initial = settings.optimality_initial;
            settingsTargetParams.optimality_final = settings.optimality_final;
            
            settingsTargetParams.epsilonInitial = obj.epsilonInitial;
            settingsTargetParams.epsilonFinal = obj.epsilonFinal;
            settingsTargetParams.epsilonPerInitial = obj.epsilonPerInitial;
            settingsTargetParams.epsilonPerFinal = obj.epsilonPerFinal;
            settingsTargetParams.epsilonIsotropyInitial = settings.epsilon_isotropy_initial;
            settingsTargetParams.epsilonIsotropyFinal = settings.epsilon_isotropy_final;
            
            
            obj.targetParamsManager = TargetParamsManager(settingsTargetParams);
            obj.targetParams = obj.targetParamsManager.targetParams;
        end
        
        function updateTargetParams(obj)
            obj.targetParamsManager.update(obj.iStep);
        end
        
        function incrementStep(obj)
            obj.iStep = obj.iStep + 1;
            if obj.shallDisplayStep
                obj.display();
            end
        end
        
        function setupEpsilons(obj,cParams,mesh)
            L = mesh.computeCharacteristicLength();
            D = mesh.computeMeanCellSize();
            obj.assignWithBackup('epsilonInitial',cParams.epsilon_initial,D);
            obj.assignWithBackup('epsilonFinal',cParams.epsilon_final,obj.epsilonInitial);
            obj.epsilonPerInitial = L;
            obj.epsilonPerFinal = obj.epsilonInitial;
        end
        
        function setWhetherShallDisplayStep(obj,settings)
            obj.shallDisplayStep = settings.printIncrementalIter;
            if isempty(settings.printIncrementalIter)
                obj.shallDisplayStep = true;
            end
        end
        
        function assignWithBackup(obj,prop,a,b)
            if ~isempty(a)
                obj.(prop) = a;
            else
                obj.(prop) = b;
            end
        end
        
    end
    
end