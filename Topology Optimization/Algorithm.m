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
            obj.optimality_tol=settings.optimality_tol;
            obj.constr_tol=settings.constr_tol;
            obj.epsilon_scalar_product_P1=settings.epsilon_scalar_product_P1;
        end 
        function sp=scalar_product(obj,f,g)
            sp=f'*(((obj.epsilon_scalar_product_P1)^2)*obj.Ksmooth+obj.Msmooth)*g;
        end
        function plotX(obj,x,physicalProblem)            
            rho_nodal=x<0;
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
    end
end
