function uBAR = DomainRedDirichledBC_2D3Dtiling(alphaBC,uBARin,NODESfaces,BasisUrb,...
    ndim,COORref,DATAONLINE,DATALIM) ;

%dbstop('4')
if nargin==0
    load('tmp.mat')
end

%%%%%%
% Compute center of each face
Nx = DATAONLINE.NdomX ;
Ny = DATAONLINE.NdomY ;

xmin = DATALIM.xmin;
xmax = DATALIM.xmax ;
ymin = DATALIM.ymin ;
ymax = DATALIM.ymax ;

if ndim == 2
    %% REference point for each subdomain  (corner, xmin, ymin)
    % ------------------------------------
    COORx = zeros(Nx*Ny,1) ; 
    COORy = zeros(Nx*Ny,1) ; 
    for ix = 1:Nx         
        for iy =1:Ny
            idom = (iy-1)*Nx + ix;
            COORx(idom) = xmin +(xmax-xmin)*(ix-1) ;
            COORy(idom) = ymin + (ymax-ymin)*(iy-1) ; 
            
        end
    end
    COOR_0 = [COORx,COORy]; 
    
    
    
    % We assume points are equally spaced on the interfaces
    % FACE x  =xminGLO
    xminTOT = xmin ;
    ymaxTOT =  ymax + (ymax-ymin)*(Ny-1) ;
    yminTOT = ymin ;
    yMID = (ymaxTOT-yminTOT)/2 ;
    xminFACE = [xminTOT,yMID]  ;
    % FACE x  =xmaxGLO
    xmaxTOT =  xmax + (xmax-xmin)*(Nx-1) ;
    xmaxFACE = [xmaxTOT,yMID]  ;
    % FACE y  =yminGLO
    xMID = (xmaxTOT-xminTOT)/2 ;
    yminFACE = [xMID,yminTOT]  ;
    % FACE y  =yminGLO
    ymaxFACE = [xMID,ymaxTOT]  ;
    REFERENCE_POINTS = {xminFACE,yminFACE,xmaxFACE,ymaxFACE}  ;
else
 DATAONLINE = DefaultField(DATAONLINE,'NdomZ',1) ; 
    Nz = DATAONLINE.NdomZ ;
    
   if Ny == 1 & Nz ==1
           REFERENCE_POINTS = [] ;

   else
    error('Option not implemented')
   end
 end

%%%%%%


uBAR = cell(size(uBARin));

for idom = 1:size(alphaBC,1)
    % Loop over domains
    FACES = find(alphaBC(idom,:)==1) ;  % FAces with imposed displacements
    
    for iface = 1:length(FACES)
        indFACE = FACES(iface) ;  % Face index
        uBAR_0 = uBARin{idom,indFACE} ;  % Imposed displacement (translation+ rotation)
        FACEnodes = NODESfaces{indFACE} ; % Nodes involved face
        COOR_FACE = COORref(FACEnodes,:) ;
         if isempty(REFERENCE_POINTS)
            
            COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
        else
            COORrefPOINT = REFERENCE_POINTS{indFACE} - COOR_0(idom,:) ;
        end
        
        COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
         BasisUrbLOC = ConstructBasisRigidBody(COORrel) ; %
        DOFface = 1:size(BasisUrbLOC,1) ;
        uBARloc = PrescribedRIGIBbodyDisp(BasisUrbLOC,DOFface,uBAR_0)  ;
        
        uBAR{idom,indFACE} = uBARloc ;
    end
    
end
%