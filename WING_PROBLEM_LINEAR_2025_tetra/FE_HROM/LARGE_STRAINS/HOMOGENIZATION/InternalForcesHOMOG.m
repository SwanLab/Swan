function  Fint = InternalForcesHOMOG(OPERFE,PoneST,DATA) ;  
% Internal forces, large strains, in terms of the Bst, W operators, as well
% as the 1st PK stress (stacked) 
% 
if nargin == 0
    load('tmp.mat')
end
nF = DATA.MESH.ndim^2 ;
for icomp = 1:nF
    icol = icomp:nF:length(PoneST) ;
    PoneST(icol) = PoneST(icol).*OPERFE.wSTs; 
end

Fint = OPERFE.BstA'*PoneST ; 

