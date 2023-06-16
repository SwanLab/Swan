function plotMesh(in)
    subMesh = in.subMesh;
    subBoundMesh = in.subBoundMesh;
    subdomains = length(subMesh);
    
    figure
    set(gcf,'units','normalized')
    set(gcf,'Position',[0, 0.1, 1, 0.9])
    
    for i=1:subdomains
        subMesh(i).plot
        subBoundMesh(i,1).mesh.plot
        if i~=subdomains
            subBoundMesh(i,2).mesh.plot
        end
    end

end