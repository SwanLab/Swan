classdef Optimizer < handle
    properties
        fhtri
        kappa
        objfunc
        stop_criteria=1;
        Msmooth
        Ksmooth
        target_parameters=struct;
        epsilon_scalar_product_P1
        shfunc_volume
        name
        niter
        optimizer
        maxiter
        physicalProblem
    end
    methods
        function obj=Optimizer(settings)
            obj.shfunc_volume=ShFunc_Volume(settings);
            obj.target_parameters=settings.target_parameters;
            obj.optimizer = settings.optimizer;
            obj.maxiter = settings.maxiter;
        end
        
        function x=solveProblem(obj,x_ini,cost,constraint,interpolation,filter)
            obj.epsilon_scalar_product_P1=1*obj.estimate_mesh_size(obj.physicalProblem.mesh.coord,obj.physicalProblem.mesh.connec);
            obj.Msmooth=obj.physicalProblem.computeMass(2);
            obj.Ksmooth=obj.physicalProblem.computeKsmooth;
            cost.computef(x_ini,obj.physicalProblem,interpolation,filter);
            constraint.computef(x_ini,obj.physicalProblem,interpolation,filter);
            obj.plotX(x_ini)
            iter=0;
            obj.print(x_ini,filter.getP0fromP1(x_ini),iter);
            while(obj.stop_criteria && iter < obj.maxiter)
                iter=iter+1;
                x=obj.updateX(x_ini,cost,constraint,interpolation,filter);
                obj.plotX(x)
                obj.print(x,filter.getP0fromP1(x),iter);
                x_ini=x;
            end
            obj.stop_criteria=1;
            obj.niter = iter;
            
        end
        
        function setPhysicalProblem(obj,pProblem)
            obj.physicalProblem = pProblem;
        end
    end
    methods (Access = private)
        function h=estimate_mesh_size(obj,coordinates,conectivities)
            x1 = coordinates(conectivities(:,1));
            x2 = coordinates(conectivities(:,2));
            x3 = coordinates(conectivities(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            h = mean(hs);
        end
        function print(obj,design_variable,design_variable_reg,iter)
            postprocess = Postprocess_TopOpt.Create(obj.optimizer);
            results.physicalVars = obj.physicalProblem.variables;
            results.design_variable = design_variable;
            results.design_variable_reg = design_variable_reg;
            postprocess.print(obj.physicalProblem,obj.name,iter,results);
        end
        function compute_physical_variables(obj)
            switch obj.physicalProblem.mesh.scale
                case 'MICRO'
                    obj.physicalProblem.computeChomog;
                case 'MACRO'
                    obj.physicalProblem.computeVariables;
            end
        end
    end
    methods (Access = protected)
        function update_physical_variables(obj,x,interpolation,filter)
            rho=filter.getP0fromP1(x);
            %Update phys problem
            matProps=interpolation.computeMatProp(rho);
            obj.physicalProblem.setMatProps(matProps);
            obj.compute_physical_variables;
        end
        function sp=scalar_product(obj,f,g)
            f = f(:);
            g = g(:);
            sp=f'*(((obj.epsilon_scalar_product_P1)^2)*obj.Ksmooth+obj.Msmooth)*g;
        end
        function plotX(obj,x)
            if any(x<0)
                rho_nodal=x<0;
            else
                rho_nodal=x;
            end
            if isempty(obj.fhtri)
                fh = figure;
                mp = get(0, 'MonitorPositions');
                select_screen=1;
                if size(mp,1) < select_screen
                    select_screen = size(mp,1);
                end
                width = mp(1,3);
                height = mp(1,4);
                size_screen_offset = round([0.7*width,0.52*height,-0.71*width,-0.611*height],0);
                set(fh,'Position',mp(select_screen,:) + size_screen_offset);
                obj.fhtri = trisurf(obj.physicalProblem.mesh.connec,obj.physicalProblem.mesh.coord(:,1),obj.physicalProblem.mesh.coord(:,2),double(rho_nodal), ...
                    'EdgeColor','none','LineStyle','none','FaceLighting','phong');
                view([0,90]);
                colormap(flipud(gray));
                set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
                drawnow;
            else
                set(obj.fhtri,'FaceVertexCData',double(rho_nodal));
                set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
                drawnow;
            end
        end
    end
    methods (Static)
        %         function physicalProblem=updateEquilibrium(x,physicalProblem,interpolation,filter)
        %             rho=filter.getP0fromP1(x);
        %             %Update phys problem
        %             matProps=interpolation.computeMatProp(rho);
        %             physicalProblem.setMatProps(matProps);
        %             physicalProblem.computeVariables;
        %         end
    end
end
