function inter=CheckPointsInsideTetrahedron(V,F_Local_new, Vertex_id)
for j=1:size(F_Local_new,1)
    P1=V(F_Local_new(j,1),:); 
    P2=V(F_Local_new(j,2),:); 
    P3=V(F_Local_new(j,3),:); 
    P4=V(Vertex_id,:);
    
    % Box around tetrahedron
    mi=[min([P1(1) P2(1) P3(1) P4(1)]) min([P1(2) P2(2) P3(2) P4(2)]) min([P1(3) P2(3) P3(3) P4(3)])];
    ma=[max([P1(1) P2(1) P3(1) P4(1)]) max([P1(2) P2(2) P3(2) P4(2)]) max([P1(3) P2(3) P3(3) P4(3)])];
    for i=1:size(V,1)
        check=(i==F_Local_new(j,1))||(i==F_Local_new(j,2))||(i==F_Local_new(j,3))||(i==Vertex_id);
        if(~check)
            inter=CheckInsideTetrahedron(P1,P2,P3,P4,V(i,1),V(i,2),V(i,3),mi,ma);
            if(inter),break; end
        end
    end
    if(inter),break; end
end
