function plotStresses(subBoundMesh,Stress,React)

    subdomains = size(subBoundMesh,1);
    for i=subdomains:-1:1

        Stress{i}.plot
        X = subBoundMesh(i,1).mesh.coord(1,1);
        title("Stresses (Boundary "+(i-1)+"-"+(i)+")"+" [X="+X+"]")
        xlabel('Stress magnitude')
        ylabel('Section position')
        legend('Normal stress $\sigma_x$','Shear stress $\tau_{xy}$','interpreter','latex')
    end
end