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
    end
    methods
        function obj=Optimizer(settings)       
            obj.shfunc_volume=ShFunc_Volume(settings);  
            obj.target_parameters=settings.target_parameters;
            obj.optimizer = settings.optimizer;
        end 
        
        function x=solveProblem(obj,x_ini,cost,constraint,physProblem,interpolation,filter) 
            obj.epsilon_scalar_product_P1=1*obj.estimate_mesh_size(physProblem.mesh.coord,physProblem.mesh.connec);
            obj.Msmooth=physProblem.computeMass(2);
            obj.Ksmooth=physProblem.computeKsmooth;
            cost.computef(x_ini,physProblem,interpolation,filter);
            constraint.computef(x_ini,physProblem,interpolation,filter);  
            obj.plotX(x_ini,physProblem)
            iter=0;
            obj.print(x_ini,physProblem,filter.getP0fromP1(x_ini),iter);
            while(obj.stop_criteria)
                iter=iter+1;
                x=obj.updateX(x_ini,cost,constraint,physProblem,interpolation,filter);
                obj.plotX(x,physProblem)
               % obj.print(x,physProblem,filter.getP0fromP1(x),iter);
                x_ini=x;                
            end
            obj.stop_criteria=1;
            obj.niter = iter;

        end
        
        function sp=scalar_product(obj,f,g)
            sp=f'*(((obj.epsilon_scalar_product_P1)^2)*obj.Ksmooth+obj.Msmooth)*g;
        end
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
        function plotX(obj,x,physicalProblem)
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
                obj.fhtri = trisurf(physicalProblem.mesh.connec,physicalProblem.mesh.coord(:,1),physicalProblem.mesh.coord(:,2),double(rho_nodal), ...
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
       
        function print(obj,design_variable,physicalProblem,design_variable_reg,iter)
            postprocess = Postprocess_TopOpt.Create(obj.optimizer);
            results.physicalVars = physicalProblem.variables;
            results.design_variable = design_variable;
            results.design_variable_reg = design_variable_reg;
            postprocess.print(physicalProblem,obj.name,iter,results);

        end
        
    end
    methods (Static)
        function physicalProblem=updateEquilibrium(x,physicalProblem,interpolation,filter)
            rho=filter.getP0fromP1(x);
            %Update phys problem
            matProps=interpolation.computeMatProp(rho);
            physicalProblem.setMatProps(matProps);
            physicalProblem.computeVariables;
        end
    end
end
