function intersect=CheckInputMesh(V,F,T)
if(exist('T','var'))
    F1=T(:,[1 2 3]);
    F2=T(:,[1 4 2]); 
    F3=T(:,[1 3 4]); 
    F4=T(:,[3 2 4]); 
    F5=F;
    F=[F1;F2;F3;F4;F5];
end

[intersect1 FV]=CheckMeshInterSections(V,F);
    
if(intersect1)
    warning('CheckInputMesh:Input','Mesh intersections detected');
end
intersect2=CheckFaceOrientations(F);
if(intersect2)
    warning('CheckInputMesh:Input','Non consistent face orientation');
end

F2=F(:);
Find=unique(F2);
intersect3=false;
for i=1:length(Find)
    if(sum(F2==Find(i))==1), intersect3=true; break;
    end
end
if(intersect3)
    warning('CheckInputMesh:Input','Mesh contains unique vertices!');
end

volume=CheckVolumeFaceMesh(V,F);
if(volume<0), intersect4=true; else intersect4=false; end


intersect=[intersect1 intersect2 intersect3 intersect4];