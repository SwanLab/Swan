classdef IncrementalScheme < handle
    
    properties (Access = public)
        settings
        incropt
        coord
        connec
        epsilon
        epsilon_initial
        epsilon0
        epsilon_isotropy
    end
    
    methods (Access = public)
        
        function obj = IncrementalScheme(settings,mesh)
            obj.settings = settings;
            nsteps = settings.nsteps;
            obj.coord = mesh.coord;
            obj.connec = mesh.connec;
            if isempty(settings.epsilon_initial)
                obj.epsilon_initial = mesh.computeMeanCellSize();
            else
                obj.epsilon_initial = settings.epsilon_initial;
            end
            obj.epsilon = obj.epsilon_initial;
            obj.epsilon0 = mesh.computeCharacteristicLength();
            obj.incropt.alpha_vol = obj.generate_incr_sequence(1/nsteps,1,nsteps,'linear');
            obj.incropt.alpha_constr = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_optimality = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon_vel = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon_per = obj.generate_incr_sequence(-1,0,nsteps,'logarithmic');
            if strcmp(obj.settings.ptype,'MICRO')
                obj.incropt.alpha_epsilon_isotropy = obj.generate_incr_sequence(0,1,nsteps,'linear');
            end
        end
        
        function update_target_parameters(obj,t,cost,constraint,optimizer)
            target_parameters.Vfrac = (1-obj.incropt.alpha_vol(t))*obj.settings.Vfrac_initial+obj.incropt.alpha_vol(t)*obj.settings.Vfrac_final;
            target_parameters.epsilon_perimeter = (1-obj.incropt.alpha_epsilon_per(t))*obj.epsilon0+obj.incropt.alpha_epsilon_per(t)*obj.epsilon;
            target_parameters.epsilon = (1-obj.incropt.alpha_epsilon(t))*obj.epsilon_initial+obj.incropt.alpha_epsilon(t)*obj.epsilon;
            target_parameters.epsilon_velocity = (1-obj.incropt.alpha_epsilon_vel(t))*obj.epsilon0+obj.incropt.alpha_epsilon_vel(t)*obj.epsilon;
            target_parameters.constr_tol = (1-obj.incropt.alpha_constr(t))*obj.settings.constr_initial+obj.incropt.alpha_constr(t)*obj.settings.constr_final;
            target_parameters.optimality_tol = (1-obj.incropt.alpha_optimality(t))*obj.settings.optimality_initial+obj.incropt.alpha_optimality(t)*obj.settings.optimality_final;
            
            if strcmp(obj.settings.ptype,'MICRO')
                target_parameters.epsilon_isotropy = (1-obj.incropt.alpha_epsilon_isotropy(t))*obj.settings.epsilon_isotropy_initial+obj.incropt.alpha_epsilon_isotropy(t)*obj.settings.epsilon_isotropy_final;
            end
            
            cost.target_parameters = target_parameters;
            constraint.target_parameters = target_parameters;
            optimizer.target_parameters = target_parameters;
        end
        
    end
    
    methods (Access = private)
        
        function x = generate_incr_sequence (obj,x1,x2,nsteps,type,factor)
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