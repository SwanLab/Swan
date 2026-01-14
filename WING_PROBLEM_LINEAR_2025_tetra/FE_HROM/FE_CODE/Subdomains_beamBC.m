function [GammaFACES,TH,uBAR,alphaBC,betaBC] = Subdomains_beamBC(Nx,Ny,FIXED_FACE,uBAR_fixed) ;  
 if nargin == 0
    load('tmp5.mat')
    Nx = 4; Ny = 3 ; 
end

%% GAMMAface, see MultiLearn.pdf and Implementation.pdf (f1=xmax, f2=ymin, f3=xmin, f4 = ymax, f5 = zmax, f6 = zmin)
% Plate geometry of dimension NDOMdir = [Nx  Ny 1]
 NDOM = Nx*Ny  ;

G = sparse(NDOM,NDOM) ; % For lines
Gpoints = sparse(NDOM,NDOM) ; % For points

% xmin -> face 1
% ymin --> face 2
% xmax --> face 3
% ymax ---> face 4

%% Connections along x directions
for iz = 1:Ny
    COLX1 = (iz-1)*Nx+(1:Nx-1);
    COLX2 = COLX1 + 1;
    IND = sub2ind(size(G),COLX1,COLX2) ;
    G(IND) = 3 ;
    IND = sub2ind(size(G),COLX2,COLX1) ;
    G(IND) = 1 ;     
end
%% Connections along z directions

for ix = 1:Nx
    COLZ1 = (ix-1)+(1:Nx:(NDOM-Nx));
    COLZ2 = COLZ1 + Nx;
    IND = sub2ind(size(G),COLZ1,COLZ2) ;
    G(IND) = 4 ;
    IND = sub2ind(size(G),COLZ2,COLZ1) ;
    G(IND) = 2 ;
end

% --- POINTS ---------------------
for iz = 1:Ny-1
    COLX1 = (iz-1)*Nx+(1:Nx-1);
    COLX2 = iz*Nx+(2:Nx);   
    IND = sub2ind(size(Gpoints),COLX1,COLX2) ;
    Gpoints(IND) = 4 ;
    IND = sub2ind(size(Gpoints),COLX2,COLX1) ;
    Gpoints(IND) = 1 ;
    
    COLX1 = (iz-1)*Nx+(2:Nx);
    COLX2 = iz*Nx+(1:Nx-1);   
    IND = sub2ind(size(Gpoints),COLX1,COLX2) ;
    Gpoints(IND) = 3 ;
    IND = sub2ind(size(Gpoints),COLX2,COLX1) ;
    Gpoints(IND) = 2 ;
end




%%% Interface faces  (TH = NDOM x nFACES)
% --------------------------------------
nFACES = 4 ;
TH = zeros(NDOM,nFACES);
for idom = 1:NDOM
    Ny = find(G(idom,:)>0) ;
    TH(idom,G(idom,Ny)) =1 ;
end

%% Dirichlet Boundary faces (along with conditions)
% -------------------------------------------------
uBAR = cell(NDOM,nFACES);
for ibnd = 1:length(FIXED_FACE)
    switch FIXED_FACE{ibnd}
        case 'xMIN'            
            COL =  (1:Nx:(NDOM));
            indFACE = 1; 
        case 'xMAX'
            ix = Nx ;
            COL = (ix-1)+(1:Nx:(NDOM));
            indFACE = 3; 
        case 'yMIN'
            COL =  (1:Nx);
            indFACE = 2 ; 
        case 'yMAX'
            iz = Ny ;
            COL = (iz-1)*Nx+(1:Nx);          
            indFACE = 4 ; 
    end    
      uBAR(COL,indFACE) = {uBAR_fixed{ibnd}} ; ;
end

[alphaBC_comple ] = cellfun(@isempty,uBAR) ; 
alphaBC = (alphaBC_comple==0) ; %

%% Free surfaces
TH_compl = (TH==0); 
betaBC = TH_compl.*alphaBC_comple ; 

GammaFACES.FACES = G ; 
GammaFACES.LINES = Gpoints; 
GammaFACES.POINTS = sparse(NDOM,NDOM) ; 