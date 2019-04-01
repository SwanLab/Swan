classdef IncrementalScheme < handle
    
    properties (Access = public)        
        epsilonInitial
        epsilonFinal
        epsilonPerInitial
        epsilonPerFinal
        epsilonVelInitial
        epsilonVelFinal
        %         epsilon_isotropy
    end
    
    properties (GetAccess = public, SetAccess = private)
        iStep
        nSteps
        targetParams = TargetParameters
    end
    
    properties (Access = private)
        targetParamsManager
        
        settings
        
        cost
        constraint
        optimizer
        
        scale
        
        shallDisplayStep
    end
    
    methods (Access = public)
        
        function obj = IncrementalScheme(settings,mesh)
            obj.init(settings,mesh);
            obj.createTargetParams(settings);
            obj.initTargetParams();
        end
        
        function link(obj,cost,constraint,optimizer)
            obj.cost = cost;
            obj.constraint = constraint;
            obj.optimizer = optimizer;
            obj.assignTargetParams();
        end
        
        function next(obj)
            obj.incrementStep();
            obj.updateTargetParams();
        end

        
        function display(obj)
            disp(['Incremental step: ',int2str(obj.iStep),' of ',int2str(obj.nSteps)]);
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
            obj.settings = settings;
            obj.iStep = 0;
            obj.nSteps = settings.nsteps;
            obj.scale = settings.ptype;
            obj.setupEpsilons(settings.epsilon_initial,mesh);
            obj.setWhetherShallDisplayStep(settings);
        end
        
        function createTargetParams(obj,settings)
            settingsTargetParams = struct;
            settingsTargetParams.nSteps = obj.nSteps;
            settingsTargetParams.scale = obj.scale;
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
            settingsTargetParams.epsilonVelInitial = obj.epsilonVelInitial;
            settingsTargetParams.epsilonVelFinal = obj.epsilonVelFinal;
            settingsTargetParams.epsilonIsotropyInitial = settings.epsilon_isotropy_initial;
            settingsTargetParams.epsilonIsotropyFinal = settings.epsilon_isotropy_final;
            
            
            obj.targetParamsManager = TargetParamsManager(settingsTargetParams);
        end
        
        function incrementStep(obj)
            obj.iStep = obj.iStep + 1;
            if obj.shallDisplayStep
                obj.display();
            end
        end
        
        function initTargetParams(obj)
            obj.computeTargetParams(1)
        end
        
        function updateTargetParams(obj)
            obj.computeTargetParams(obj.iStep)
        end
        
        function computeTargetParams(obj,iStep)
            obj.targetParamsManager.update(iStep);
            
            obj.targetParams.Vfrac = obj.targetParamsManager.volumeFrac.value;
            obj.targetParams.epsilon = obj.targetParamsManager.epsilon.value;
            obj.targetParams.epsilon_velocity = obj.targetParamsManager.epsilonVel.value;
            obj.targetParams.epsilon_perimeter = obj.targetParamsManager.epsilonPer.value;
            obj.targetParams.constr_tol = obj.targetParamsManager.constraintTol.value;
            obj.targetParams.optimality_tol = obj.targetParamsManager.optimalityTol.value;
            if strcmp(obj.settings.ptype,'MICRO')
                obj.targetParams.epsilon_isotropy = obj.targetParamsManager.epsilonIsotropy.value;
            end
        end
        
        function assignTargetParams(obj)
            obj.cost.target_parameters = obj.targetParams;
            obj.constraint.target_parameters = obj.targetParams;
            obj.optimizer.target_parameters = obj.targetParams;
        end
        
        function setupEpsilons(obj,initialEpsilon,mesh)
            if ~isempty(initialEpsilon)
                obj.epsilonInitial = initialEpsilon;
            else
                obj.epsilonInitial = mesh.computeMeanCellSize();
            end
            obj.epsilonFinal = obj.epsilonInitial;
            obj.epsilonPerInitial = mesh.computeCharacteristicLength();
            obj.epsilonVelInitial = mesh.computeCharacteristicLength();
            obj.epsilonPerFinal = obj.epsilonInitial;
            obj.epsilonVelFinal = obj.epsilonInitial;
        end
        
        function itShall = setWhetherShallDisplayStep(obj,settings)
            itShall = settings.printIncrementalIter;
            if isempty(settings.printIncrementalIter)
                itShall = true;
            end
        end
        
    end
    
end