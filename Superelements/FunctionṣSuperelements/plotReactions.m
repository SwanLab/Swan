function plotReactions(reactions,subBoundMesh,sp)
        X = subBoundMesh.mesh.coord(1,1);
        Y = subBoundMesh.mesh.coord(:,2);
        
        figure
        plot(reactions(:,1),Y,'.b')
        hold on
        plot([0 0],[0 0.5],'-')
        quiver(zeros(size(Y(sp:sp:end-sp))),Y(sp:sp:end-sp),reactions(sp:sp:end-sp,1),zeros(size(reactions(sp:sp:end-sp,1))),"off",'b')
        title("Reaction X [X="+X+"]")
        xlabel('Force magnitude')
        ylabel('Section position')

        figure
        plot(reactions(:,2),Y,'.b')
        hold on
        plot([0 0],[0 0.5],'-')
        quiver(zeros(size(Y(sp:sp:end-sp))),Y(sp:sp:end-sp),reactions(sp:sp:end-sp,2),zeros(size(reactions(sp:sp:end-sp,2))),"off",'b')
        title("Reaction Y [X="+X+"]")
        xlabel('Force magnitude')
        ylabel('Section position')
end