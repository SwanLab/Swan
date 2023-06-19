function plotMesh(in)
    subMesh = in.subMesh;
    subBoundMesh = in.subBoundMesh;
    subdomains = length(subMesh);
    
    figure
    set(gcf,'units','normalized')
    title('Mesh decomposed')
    xlabel('Length (L)')
    ylabel('Height (H)')
    
    for i=1:subdomains
        subMesh(i).plot
        subBoundMesh(i,1).mesh.plot
        if i~=subdomains
            subBoundMesh(i,2).mesh.plot
        end
    end

end