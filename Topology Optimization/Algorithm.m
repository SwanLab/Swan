classdef Algorithm < handle
    properties
        fhtri
        shfunc_volume
        Msmooth
        Ksmooth
        constr_tol
        optimality_tol
        epsilon_scalar_product_P1;
    end
    methods
        function obj=Algorithm(settings)
            obj.shfunc_volume=ShFunc_Volume(settings.volume);
            obj.epsilon_scalar_product_P1=settings.epsilon_scalar_product_P1;
        end
        function sp=scalar_product(obj,f,g)
            sp=f'*(((obj.epsilon_scalar_product_P1)^2)*obj.Ksmooth+obj.Msmooth)*g;
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
            else
                set(obj.fhtri,'FaceVertexCData',double(rho_nodal));
                set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
                drawnow;
            end
        end     
        function printX(obj,name,x,physicalProblem,iter)
            postprocess = Postprocess;
            postprocess.ToGiD(name,physicalProblem,iter);
            postprocess.ToGiDpostX(name,x,physicalProblem,iter);
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
