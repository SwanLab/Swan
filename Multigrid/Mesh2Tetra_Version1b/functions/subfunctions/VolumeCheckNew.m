function [inter,F2,nF2,Tnew]=VolumeCheckNew(V,F,nF,T,nT,Frows,F_Local,F_Local_new, Vertex_id,Vertex_id_Local,Volume_Original)
nTold=nT;
[V,F,nF,T,nT]=process(V,F,nF,T,nT,Frows,F_Local,F_Local_new, Vertex_id);
F2=F; nF2=nF; 
inter=VolumeCheck(V,F,nF,T,nT,Volume_Original);
Tnew=T(nTold+1:nT,:);
