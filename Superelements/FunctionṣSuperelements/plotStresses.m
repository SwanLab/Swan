function plotStresses(subBoundMesh,stress,type)
    if type == "Two"
        for i=1:length(stress)
            stress{i}.plot
            X = subBoundMesh(i,1).mesh.coord(1,1);
            title("Stresses" + " [X="+X+"]")
            xlabel('Stress magnitude')
            ylabel('Section position')
            legend('Normal stress $\sigma_x$','Shear stress $\tau_{xy}$','interpreter','latex')
        end
    elseif type == "Three" 
        for i=1:length(stress)
            stress{i}.plot
            if rem(i,2)==0
               X = subBoundMesh(round(i/2),2).mesh.coord(1,1);
            else
               X = subBoundMesh(round(i/2),1).mesh.coord(1,1);
            end

            title("Stresses" + " [X="+X+"]")
            xlabel('Stress magnitude')
            ylabel('Section position')
            legend('Normal stress $\sigma_x$','Shear stress $\tau_{xy}$','interpreter','latex')
        end
    end

end