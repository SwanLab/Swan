classdef Optimizer < handle
    %Optimizer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stop_updating = 0;
        stop_vars
        target_parameters = struct;
        nconstr
        constraint_case
        ini_design_value = 1;
        hole_value = 0;
        holes
        postprocess
    end
    
    properties (Access = ?Optimizer_Constrained)
        case_file
        printing
        monitoring
        mesh
    end
    
    methods
        function obj = Optimizer(settings)
            obj.nconstr = settings.nconstr;
            obj.case_file=settings.case_file;
            obj.holes.N_holes = settings.N_holes;
            obj.holes.R_holes = settings.R_holes;
            obj.holes.phase_holes = settings.phase_holes;
            obj.target_parameters = settings.target_parameters;
            obj.constraint_case=settings.constraint_case;
            obj.postprocess = Postprocess_TopOpt.Create(settings.optimizer);
        end
        
        function x = compute_initial_design(obj,initial_case,optimizer)
            geometry = Geometry(obj.mesh,'LINEAR');
            x = obj.ini_design_value*ones(geometry.interpolation.npnod,1);
            switch initial_case
                case 'circle'
                    W = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
                    H = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
                    radius = 0.2*min([W,H]);
                    initial_holes = (obj.mesh.coord(:,1)-center_x).^2 + (obj.mesh.coord(:,2)-center_y).^2 - radius^2 < 0;
                    x(initial_holes) = obj.hole_value;
                    
                case 'horizontal'
                    initial_holes = obj.mesh.coord(:,2) > 0.6 | obj.mesh.coord(:,2) < 0.4;
                    x(initial_holes) = obj.hole_value;
                    
                case 'square'
                    W = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
                    H = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.mesh.coord(:,1)) + min(obj.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.mesh.coord(:,2)) + min(obj.mesh.coord(:,2)));
                    offset_x = 0.2*W;
                    offset_y = 0.2*H;
                    xrange = obj.mesh.coord(:,1) < (center_x+offset_x) & obj.mesh.coord(:,1) > (center_x-offset_x);
                    yrange = obj.mesh.coord(:,2) < (center_y+offset_y) & obj.mesh.coord(:,2) > (center_y-offset_y);
                    initial_holes = and(xrange,yrange);
                    x(initial_holes) = obj.hole_value;
                    
                case 'feasible'
                    initial_holes = false(size(obj.mesh.coord,1),1);
                    x(initial_holes) = obj.hole_value;
                    
                case 'rand'
                    initial_holes = rand(size(obj.mesh.coord,1),1) > 0.1;
                    x(initial_holes) = obj.hole_value;
                    
                case 'holes'
                    L(1) = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
                    L(2) = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));
                    switch obj.mesh.pdim
                        case '2D'
                            initial_holes = ceil(max(cos((obj.holes.N_holes(2)+1)*(obj.mesh.coord(:,2)*pi)/L(2)+obj.holes.phase_holes(2)).*cos((obj.holes.N_holes(1)+1)*(obj.mesh.coord(:,1)*pi)/L(1)+obj.holes.phase_holes(1)) + obj.holes.R_holes-1 ,0))>0;
                        case '3D'
                            L(3) = max(obj.mesh.coord(:,3)) - min(obj.mesh.coord(:,3));
                            initial_holes = ceil(max(cos((obj.holes.N_holes(3)+1)*(obj.mesh.coord(:,3)*pi)/L(3)+obj.holes.phase_holes(3)).*cos((obj.holes.N_holes(2)+1)*(obj.mesh.coord(:,2)*pi)/L(2)+obj.holes.phase_holes(2)).*cos((obj.holes.N_holes(1)+1)*(obj.mesh.coord(:,1)*pi)/L(1)+obj.holes.phase_holes(1)) + obj.holes.R_holes-1 ,0))>0;
                    end
                    x(initial_holes) = obj.hole_value;
                    
                    bc = unique([obj.mesh.dirichlet(:,1); obj.mesh.pointload(:,1)]);
                    if any(x(bc)>0)
                        warning('At least one BC is set on a hole')
                    end
                case 'full'
                otherwise
                    error('Invalid initial value of design variable.');
            end
            %% !! PROVISIONAL !!
            if strcmp(optimizer,'SLERP') %|| strcmp(optimizer,'HAMILTON-JACOBI')
                sqrt_norma = obj.optimizer_unconstr.scalar_product.computeSP(x,x);
                x = x/sqrt(sqrt_norma);
            end
        end
    end
    
    methods (Abstract)
        % x = updateX(obj,x_ini,cost,constraint); %% !! IPOPT doesn't use it (black box) !!
    end
    
    methods (Access = protected)
        function cons = setConstraint_case(obj,constraint)
            cons = constraint;
            switch obj.constraint_case
                case 'EQUALITY'
                case 'INEQUALITY'
                    contr_inactive_value = -obj.objfunc.lambda(:)./obj.objfunc.penalty(:);
                    inactive_constr = contr_inactive_value' > constraint.value;
                    cons.value(inactive_constr) = contr_inactive_value(inactive_constr);
                    cons.gradient(:,inactive_constr) = 0;
                otherwise
                    error('Constraint case not valid.');
            end
        end
        function print(obj,design_variable,iter)
            if ~(obj.printing)
                return
            end
            
            %results.physicalVars = obj.physicalProblem.variables;
            results.design_variable = design_variable;
            results.iter = iter;
            results.case_file = obj.case_file;
            %results.design_variable_reg = design_variable_reg;
            obj.postprocess.print(obj.mesh,results);
        end
        function writeToFile(obj,nstep,cost,constraint)
            if ~(obj.printing)
                return
            end
            if obj.niter==1
                msh_file = fullfile('Output',obj.case_file,strcat(obj.case_file,'.txt'));
                fid_mesh = fopen(msh_file,'wt');
            else
                msh_file = fullfile('Output',obj.case_file,strcat(obj.case_file,'.txt'));
                fid_mesh = fopen(msh_file,'at');
            end
            fprintf(fid_mesh,'-----------------------------------------------------------------------------------------------\n');
            fprintf(fid_mesh,'\n');
            fprintf(fid_mesh,'Iteration: %i \n',obj.niter);
            fprintf(fid_mesh,'Nstep: %i \n',nstep);
            fprintf(fid_mesh,'Cost  %f \n',cost.value);
            for i=1:length(cost.ShapeFuncs)
                fprintf(fid_mesh,strcat('-Cost ',num2str(i),': %f \n'),cost.ShapeFuncs{i}.value);
            end
            fprintf(fid_mesh,'Constraint: %f \n',constraint.value);
            for i=1:length(constraint.ShapeFuncs)
                fprintf(fid_mesh,strcat('-Constraint ',num2str(i),': %f \n'),constraint.ShapeFuncs{i}.value);
            end
            
            switch obj.optimizer
                case {'SLERP','PROJECTED GRADIENT','HAMILTON-JACOBI'}
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.optimizer_unconstr.opt_cond);
                    fprintf(fid_mesh,'Kappa: %f \n',obj.optimizer_unconstr.kappa);
                case 'MMA'
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.kktnorm);
                case 'IPOPT'
                    fprintf(fid_mesh,'Optimality tolerance: %f \n',obj.data.inf_du);
            end
            fprintf(fid_mesh,'\n');
            fclose(fid_mesh);
        end
    end
end

