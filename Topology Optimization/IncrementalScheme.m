classdef IncrementalScheme < handle
    
    properties (Access = public)
        incropt
        
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
        targetParams
    end
    
    properties (Access = private)
        minEpsilon
        
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
            obj.generateIncrementalSequences();
            obj.initTargetParams();
        end
        
        function link(obj,cost,constraint,optimizer)
            obj.cost = cost;
            obj.constraint = constraint;
            obj.optimizer = optimizer;
        end
        
        function next(obj)
            obj.incrementStep();
            obj.updateTargetParams();
            obj.assignTargetParams();
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
            obj.computeMinimumEpsilon(mesh);
            obj.setupEpsilons(settings.epsilon_initial);
            obj.setWhetherShallDisplayStep(settings);
        end
        
        function generateIncrementalSequences(obj)
            nsteps = obj.nSteps;
            
            obj.incropt.volumeFrac = LinearSequence(1/nsteps,1,nsteps,obj.settings.Vfrac_initial,obj.settings.Vfrac_final);
            obj.incropt.constraintTol = LinearSequence(0,1,nsteps,obj.settings.constr_initial,obj.settings.constr_final);
            obj.incropt.optimalityTol = LinearSequence(0,1,nsteps,obj.settings.optimality_initial,obj.settings.optimality_final);
            obj.incropt.epsilon = LinearSequence(0,1,nsteps,obj.epsilonInitial,obj.epsilonFinal);
            obj.incropt.epsilonVel = LinearSequence(0,1,nsteps,obj.epsilonVelInitial,obj.epsilonVelFinal);
            obj.incropt.epsilonPer = LogarithmicSequence(-1,0,nsteps,obj.epsilonPerInitial,obj.epsilonPerFinal);
            
            %             obj.incropt.epsilonPer = EpsilonSequence(1/nsteps,1,nsteps,obj.epsilonPerInitial,obj.epsilonPerFinal);
            %             obj.incropt.epsilonPer = CustomSequence(0,1,nsteps,obj.epsilonPerInitial,obj.epsilonPerFinal);
            obj.incropt.epsilonPer = FreeSequence(0,1,nsteps,obj.epsilonPerInitial,obj.epsilonPerFinal);
            
            if strcmp(obj.scale,'MICRO')
                obj.incropt.epsilonIsotropy = LinearSequence(0,1,nsteps,obj.settings.epsilon_isotropy_initial,obj.settings.epsilon_isotropy_final);
            end
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
            obj.incropt.volumeFrac.update(iStep);
            obj.incropt.constraintTol.update(iStep);
            obj.incropt.optimalityTol.update(iStep);
            obj.incropt.epsilon.update(iStep);
            obj.incropt.epsilonVel.update(iStep);
            obj.incropt.epsilonPer.update(iStep);
            
            obj.targetParams.Vfrac = obj.incropt.volumeFrac.value;
            obj.targetParams.epsilon = obj.incropt.epsilon.value;
            obj.targetParams.epsilon_velocity = obj.incropt.epsilonVel.value;
            obj.targetParams.epsilon_perimeter = obj.incropt.epsilonPer.value;
            obj.targetParams.constr_tol = obj.incropt.constraintTol.value;
            obj.targetParams.optimality_tol = obj.incropt.optimalityTol.value;
            if strcmp(obj.settings.ptype,'MICRO')
                obj.targetParams.epsilon_isotropy = obj.incropt.epsilonIsotropy.value;
            end
        end
        
        function assignTargetParams(obj)
            obj.cost.target_parameters = obj.targetParams;
            obj.constraint.target_parameters = obj.targetParams;
            obj.optimizer.target_parameters = obj.targetParams;
        end
        
        function setupEpsilons(obj,initialEpsilon)
            if ~isempty(initialEpsilon)
                obj.epsilonInitial = initialEpsilon;
            else
                obj.epsilonInitial = obj.minEpsilon;
            end
            obj.epsilonFinal = obj.epsilonInitial;
            obj.epsilonPerInitial = obj.minEpsilon;
            obj.epsilonVelInitial = obj.minEpsilon;
            obj.epsilonPerFinal = obj.epsilonFinal;
            obj.epsilonVelFinal = obj.epsilonFinal;
        end
        
        function computeMinimumEpsilon(obj,mesh)
            obj.minEpsilon = mesh.computeCharacteristicLength();
        end
        
        function itShall = setWhetherShallDisplayStep(obj,settings)
            itShall = settings.printIncrementalIter;
            if isempty(settings.printIncrementalIter)
                itShall = true;
            end
        end
        
    end
    
end