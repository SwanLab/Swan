classdef Incremental_Scheme < handle
    properties
        settings
        incropt
        coord
        connec
        epsilon
        epsilon_initial
        epsilon0
        epsilon_isotropy
    end
    methods
        function obj=Incremental_Scheme(settings, mesh)
            obj.settings=settings;
            nsteps=settings.nsteps;
            obj.coord=mesh.coord;
            obj.connec=mesh.connec;
            obj.epsilon=obj.estimate_mesh_size;
            if isempty(settings.epsilon_initial)
                obj.epsilon_initial=obj.epsilon;
            end
            obj.incropt.alpha_vol = obj.generate_incr_sequence(1/nsteps,1,nsteps,'linear');
            obj.incropt.alpha_constr = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_optimality= obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon=obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_epsilon_per = obj.generate_incr_sequence(-1,0,nsteps,'logarithmic');
            if strcmp(obj.settings.ptype,'MICRO')
                obj.incropt.alpha_epsilon_isotropy=obj.generate_incr_sequence(0,1,nsteps,'linear');
            end
        end
        function update_target_parameters(obj,t,cost, constraint, optimizer)
            target_parameters.Vfrac = (1-obj.incropt.alpha_vol(t))*obj.settings.Vfrac_initial+obj.incropt.alpha_vol(t)*obj.settings.Vfrac_final;
            target_parameters.epsilon_perimeter = (1-obj.incropt.alpha_epsilon_per(t))*obj.epsilon0+obj.incropt.alpha_epsilon_per(t)*obj.epsilon;
            target_parameters.epsilon = (1-obj.incropt.alpha_epsilon(t))*obj.epsilon_initial+obj.incropt.alpha_epsilon(t)*obj.epsilon;
            target_parameters.constr_tol = (1-obj.incropt.alpha_constr(t))*obj.settings.constr_initial+obj.incropt.alpha_constr(t)*obj.settings.constr_final;
            target_parameters.optimality_tol = (1-obj.incropt.alpha_optimality(t))*obj.settings.optimality_initial+obj.incropt.alpha_optimality(t)*obj.settings.optimality_final;
           
            % TO BE DISCUSSED IN MEETING %%%%%%%%
            if strcmp(obj.settings.ptype,'MICRO')
                target_parameters.epsilon_isotropy = (1-obj.incropt.alpha_epsilon_isotropy(t))*obj.settings.epsilon_isotropy_initial+obj.incropt.alpha_epsilon_isotropy(t)*obj.settings.epsilon_isotropy_final;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            cost.target_parameters=target_parameters;
            constraint.target_parameters=target_parameters;
            optimizer.target_parameters=target_parameters;
        end
        function h=estimate_mesh_size(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            obj.epsilon0 = norm(xmax-xmin)/2;
            
            x1 = obj.coord(obj.connec(:,1));
            x2 = obj.coord(obj.connec(:,2));
            x3 = obj.coord(obj.connec(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            h = mean(hs);
        end
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