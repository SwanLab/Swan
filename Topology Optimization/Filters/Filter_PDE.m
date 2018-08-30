classdef Filter_PDE < Filter
    properties
        dvolu
        rhs
        A_nodal_2_gauss
    end
    
    methods
        function obj = Filter_PDE(problemID,scale)
            obj@Filter(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj);
            obj.P_operator=obj.computePoperator(obj.diffReacProb.element.M);
            obj.dvolu = sparse(1:obj.diffReacProb.geometry.interpolation.nelem,1:obj.diffReacProb.geometry.interpolation.nelem,...
                sum(obj.diffReacProb.geometry.dvolu,2));
            obj.A_nodal_2_gauss = obj.computeA;
        end
        
        function x_reg = getP1fromP1(obj,x)
            rhs_x = obj.integrate_L2_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_reg = getP1fromP0(obj,x)
            rhs_x = obj.integrate_P1_function_with_shape_function(x);
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function x_gp = getP0fromP1(obj,x)
            x_reg =  obj.getP1fromP1(x);
            x_gp = obj.A_nodal_2_gauss*x_reg;
        end
        
        % !! For SHAPE OPTIMIZATION (regularize) !!
        function x_reg = regularize(obj,x,F)
            rhs_x = obj.integrate_facet_with_shape_function(x,F);
%             x(x==0) = -1e-15;
%             S = integrate_facet_to_compute_surface(obj,x);
%             A_A0 = S/(4*pi)
%             
%             load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
%             for n = 1:length(rhs_x)
%                 b(b1(n,1),b1(n,2),b1(n,3)) = rhs_x(n);
%             end
%             
%             figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO- b')
%             subplot(2,2,1), surf(-b(:,:,2)), title('b - Root')
%             subplot(2,2,2), surf(-b(:,:,end-1)), title('b - Tip')
%             subplot(2,2,3), surf(permute(-b(ceil(size(b,1)/2),:,:),[2 3 1])), title('b - XY')
%             subplot(2,2,4), surf(permute(-b(:,ceil(size(b,2)/2),:),[1 3 2])), title('b - XZ')
            
            x_reg = obj.solve_filter(rhs_x);
        end
        
        function rhs = integrate_P1_function_with_shape_function(obj,x)
            gauss_sum=0;
            for igauss=1:size(obj.M0,2)
                gauss_sum=gauss_sum+obj.A_nodal_2_gauss'*obj.M0{igauss}*x(:,igauss);
            end
            rhs = gauss_sum;
        end
        
        function x_reg = solve_filter(obj,rhs_x)
            obj.diffReacProb.computeVariables(rhs_x);
            x_reg = obj.diffReacProb.variables.x;
        end
        
        function obj = updateEpsilon(obj,epsilon)
            obj.diffReacProb.setEpsilon(epsilon);
        end
    end
end