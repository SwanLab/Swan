function [T,Tout]=RemoveInvalidTetrahedrons(T,V,F)

NZ=zeros([size(V,1) size(V,1)],'int8');
Fdel=true([size(T,1) 1]);
for i=1:size(T,1)
    Edge=[T(i,1) T(i,2)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end

    Edge=[T(i,1) T(i,3)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end
 
    Edge=[T(i,1) T(i,4)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end
 
    Edge=[T(i,2) T(i,3)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end
 
    Edge=[T(i,2) T(i,4)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end
 
    Edge=[T(i,3) T(i,4)];
    if(NZ(Edge(1),Edge(2))==0), if(CheckVisiblePoint3D(V,F,Edge(1),Edge(2))), NZ(Edge(1),Edge(2))=1; else  NZ(Edge(1),Edge(2))=-1;  end; end
    if(NZ(Edge(1),Edge(2))<0), continue; end
    
    Fdel(i)=false;
end
Tout=T;
T(Fdel,:)=[];
Tout(~Fdel,:)=[];

