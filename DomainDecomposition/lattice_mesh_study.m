clear all
close all

filename='lattice_ex1';
a.fileName=filename;
femD = FemDataContainer(a);

boundary_mesh=femD.mesh.createBoundaryMesh();
figure(1)
for ib=1:length(boundary_mesh)
    bm=boundary_mesh{ib};
    bm.mesh.plot
    hold on
end

%% Repeat cells

nrve_dir=[3 3]; %number repetition x and y
delta= [max(femD.mesh.coord(:,1)) max(femD.mesh.coord(:,2)) ];
% dom_mesh(1,1:nrve_dir(1))=m;
% dom_mesh(2,1:nrve_dir(2))=m;
coord0  = femD.mesh.coord;
connec0 = femD.mesh.connec;
figure(2)
for jdom=1:nrve_dir(2)
    %extrude in x
    for idom=1:nrve_dir(1)
        coordIJ(:,1) = coord0(:,1)+delta(1)*(idom-1);
        coordIJ(:,2) = coord0(:,2)+delta(2)*(jdom-1);            
        newConnec    = femD.mesh.nnodes*(nrve_dir(1)*(jdom-1)+idom-1);
        connecIJ     = connec0 + newConnec;
        s.coord  = coordIJ;
        s.connec = connec0;
        m = Mesh(s);
        m.plot();
        hold on
      %  m_aux.mesh.setCoord(aux_coord)
      %  m_aux.mesh.setConnec(aux_connec)
%         m_aux.mesh.plot
      %  plot(aux_coord(:,1),aux_coord(:,2),'.')
      %  hold on
        meshes{jdom,idom}=m;
%         domain2(nrve_dir(1)*(jdom-1)+idom)=m_aux;
%         dom_mesh(jdom,idom).coord(:,1)=m_aux.mesh.coord(:,1)+delta(1)*(idom-1);
%         dom_mesh(jdom,idom).coord(:,2)=m_aux.mesh.coord(:,2)+delta(1)*(jdom-1);
%         dom_mesh(jdom,idom).connec=m_aux.mesh.connec+m_aux.mesh.nnodes*(nrve_dir(1)*(jdom-1)+idom-1);
    end
end

figure
for j=1:nrve_dir(2)
    for i=1:nrve_dir(2)  
        domain{i,j}.mesh.plot
        plot(domain{i,j}.mesh.coord(:,1),domain{i,j}.mesh.coord(:,2),'.')
        hold on
    end
end


%% rpeat cells joining domains

nrve_dir=[3 3]; %number repetition x and y
delta= [max(femD.mesh.coord(:,1)) max(femD.mesh.coord(:,2)) ];
% dom_mesh(1,1:nrve_dir(1))=m;
% dom_mesh(2,1:nrve_dir(2))=m;
coord0=femD.mesh.coord;
connec0=femD.mesh.connec;
m_aux=femD;
figure
for jdom=1:nrve_dir(2)
    %extrude in x
    for idom=1:nrve_dir(1)
        coordIJ(:,1)=coord0(:,1)+delta(1)*(idom-1);
        coordIJ(:,2)=coord0(:,2)+delta(2)*(jdom-1);
        conbecIJ=connec0+femD.mesh.nnodes*(nrve_dir(1)*(jdom-1)+idom-1);
        m_aux.mesh.setCoord(coordIJ)
        m_aux.mesh.setConnec(conbecIJ)
%         m_aux.mesh.plot
        plot(coordIJ(:,1),coordIJ(:,2),'.')
        hold on
        domain{jdom,idom}=m_aux;
%         domain2(nrve_dir(1)*(jdom-1)+idom)=m_aux;
%         dom_mesh(jdom,idom).coord(:,1)=m_aux.mesh.coord(:,1)+delta(1)*(idom-1);
%         dom_mesh(jdom,idom).coord(:,2)=m_aux.mesh.coord(:,2)+delta(1)*(jdom-1);
%         dom_mesh(jdom,idom).connec=m_aux.mesh.connec+m_aux.mesh.nnodes*(nrve_dir(1)*(jdom-1)+idom-1);
    end
end




%extrude in x
for idom=1:nrve_dir(1)
    m_aux=femD;
%     m_aux.mesh.coord(:,1)=m_aux.mesh.coord(:,1)+delta(1);
%     m_aux.mesh.conec=m_aux.mesh.connec+nnode;
%     dom_mesh(idom)=m_aux;
    dom_mesh(1,idom).coord(:,1)=m_aux.mesh.coord(:,1)+delta(1)*(idom-1);
    dom_mesh(1,idom).coord(:,2)=m_aux.mesh.coord(:,2);
    dom_mesh(1,idom).connec=m_aux.mesh.connec+m_aux.mesh.nnodes*(idom-1);
end

% extrude in y
for jdom=2:nrve_dir(2)
    m_aux=femD;
%     m_aux.mesh.coord(:,1)=m_aux.mesh.coord(:,1)+delta(1);
%     m_aux.mesh.conec=m_aux.mesh.connec+nnode;
%     dom_mesh(idom)=m_aux;
    dom_mesh(jdom,:).coord(:,1)=dom_mesh(1,:).coord(:,1);
    dom_mesh(jdom,:).coord(:,2)=dom_mesh(1,:).coord(:,2)+delta(2)*(jdom-1);
    dom_mesh(jdom,:).connec=dom_mesh(1,:).connec+m_aux.mesh.nnodes*(jdom-1);
end






