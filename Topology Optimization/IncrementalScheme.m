classdef IncrementalScheme < handle
    
    properties (Access = public)
        incropt
        
        epsilon
        epsilon_initial
        epsilon0
        epsilon_isotropy
    end
    
    properties (GetAccess = public, SetAccess = private)
        iStep
        nSteps
    end
    
    properties (Access = private)
        settings
        targetParams
        
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
            obj.setupEpsilons(settings.epsilon_initial,mesh);
            obj.setWhetherShallDisplayStep(settings);
        end
        
        function generateIncrementalSequences(obj)
            nsteps = obj.nSteps;
            
            obj.incropt.alpha_vol = obj.generateIncrementalSequence(1/nsteps,1,nsteps,'linear');
            obj.incropt.alpha_constr = obj.generateIncrementalSequence(0,1,nsteps,'linear');
            obj.incropt.alpha_optimality = obj.generateIncrementalSequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon = obj.generateIncrementalSequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon_vel = obj.generateIncrementalSequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon_per = obj.generateIncrementalSequence(-1,0,nsteps,'logarithmic');
            if strcmp(obj.scale,'MICRO')
                obj.incropt.alpha_epsilon_isotropy = obj.generateIncrementalSequence(0,1,nsteps,'linear');
            end
        end
        
        function incrementStep(obj)
            obj.iStep = obj.iStep + 1;
            if obj.shallDisplayStep
                obj.display();
            end
        end
        
        function updateTargetParams(obj)
            t = obj.iStep;
            obj.targetParams.Vfrac = (1-obj.incropt.alpha_vol(t))*obj.settings.Vfrac_initial+obj.incropt.alpha_vol(t)*obj.settings.Vfrac_final;
            obj.targetParams.epsilon_perimeter = (1-obj.incropt.alpha_epsilon_per(t))*obj.epsilon0+obj.incropt.alpha_epsilon_per(t)*obj.epsilon;
            obj.targetParams.epsilon = (1-obj.incropt.alpha_epsilon(t))*obj.epsilon_initial+obj.incropt.alpha_epsilon(t)*obj.epsilon;
            obj.targetParams.epsilon_velocity = (1-obj.incropt.alpha_epsilon_vel(t))*obj.epsilon0+obj.incropt.alpha_epsilon_vel(t)*obj.epsilon;
            obj.targetParams.constr_tol = (1-obj.incropt.alpha_constr(t))*obj.settings.constr_initial+obj.incropt.alpha_constr(t)*obj.settings.constr_final;
            obj.targetParams.optimality_tol = (1-obj.incropt.alpha_optimality(t))*obj.settings.optimality_initial+obj.incropt.alpha_optimality(t)*obj.settings.optimality_final;
            
            if strcmp(obj.settings.ptype,'MICRO')
                obj.targetParams.epsilon_isotropy = (1-obj.incropt.alpha_epsilon_isotropy(t))*obj.settings.epsilon_isotropy_initial+obj.incropt.alpha_epsilon_isotropy(t)*obj.settings.epsilon_isotropy_final;
            end
        end
        
        function assignTargetParams(obj)
            obj.cost.target_parameters = obj.targetParams;
            obj.constraint.target_parameters = obj.targetParams;
            obj.optimizer.target_parameters = obj.targetParams;
        end
        
        function itShall = setWhetherShallDisplayStep(obj,settings)
            itShall = settings.printIncrementalIter;
            if isempty(settings.printIncrementalIter)
                itShall = true;
            end
        end
        
        function setupEpsilons(obj,initialEpsilon,mesh)
            if ~isempty(initialEpsilon)
                obj.epsilon_initial = initialEpsilon;
            else
                obj.epsilon_initial = mesh.computeMeanCellSize();
            end
            obj.epsilon = obj.epsilon_initial;
            obj.epsilon0 = mesh.computeCharacteristicLength();
        end
        
        function x = generateIncrementalSequence(obj,x1,x2,nsteps,type,factor)
            switch type
                case 'linear'
                    x = linspace(x1,x2,nsteps);
                    
                case 'epsilon_sequence'
                    frac = 2;
                    kmax = ceil(log10(x1/x2)/log10(frac));
                    x = obj.epsilon0./frac.^(1:kmax);
                    
                case 'logarithmic'
                    x = logspace(x1,x2,nsteps);
                    
                case 'custom'
                    if nsteps < 2
                        x = x2;
                    else
                        isteps = 0:nsteps-1;
                        x = 1-(1-isteps/(nsteps-1)).^(factor);
                        x = (x2-x1)*x + x1;
                    end
                case 'free'
                    x = zeros(1,nsteps);
                    x(end) = 1;
                otherwise
                    error('Incremental sequence type not detected.')
            end
            
        end
        
    end
    
end