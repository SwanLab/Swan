function [stressDATA,TRANSF_COORD,stressesOUT]  = StressStrains_EIFE(COOR,CN,TypeElement,...
    d,DATA,PROPMAT,MaterialType,TRANSF_COORD,Bmat,WEIGHTSinteg)
%%%% Computation stresses EIFE method
% JAHO, 12-MARCH-2023
% ---------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
% nstrain = size(celasglo,1) ;
% TypeIntegrand = 'K';
% [weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% ngaus = length(weig) ;
% if nstrain == 3
%     nstrain=4 ;
% end
% switch typePROBLEM
%     case 'pstrain'
%         TP = 1;
%     case 'pstress'
%         TP = 2 ;
%     case '3D'
%         TP = 3 ;
% end

nelem = length(MaterialType) ;
stressAVERAGEglo = cell(nelem,1);  % Average stresses (usinc CECM points)
stressMAXglo = cell(nelem,1);  % Average stresses (usinc CECM points)

stressesREF = cell(nelem,1) ; % Non-rotated stresses
stressesROTATED = cell(nelem,1) ; 

for e = 1:nelem
    %     % Elasticity matrix of element "e"
    %     celas = celasglo(:,:,e) ;
    %     celas3Dinv = celasgloINV(:,:,e) ;
    %     celas3D = inv(celas3Dinv) ;
    % Coordinates of the nodes of element "e"
    CNlocNOD = CN(e,:) ;
    %     Xe = COOR(CNlocNOD,:)' ;
    % Displacement at   nodes of element e
    CNloc = Nod2DOF(CNlocNOD,ndim) ;
    dE = d(CNloc) ; % COARSE-SCALE DISPLACEMENTS
    strainE = Bmat{e}*dE ;
    IndexDomainLOC =  TRANSF_COORD{e}.IndexParentDomain ;
    EIFEoper = PROPMAT(MaterialType(e)).EIFE_prop(IndexDomainLOC) ;
    weig = WEIGHTSinteg.INTforces{e} ;
    ngaus = length(weig) ;
    nstrain = EIFEoper.MESH.nstrain ;
    
    stressAVERAGE = zeros(nstrain,1) ;
    stressE    = zeros(nstrain,ngaus) ;
    stressesREF{e} = zeros(nstrain,ngaus) ; 
    for  g = 1:ngaus
        istrain = small2large(g,nstrain ) ;
        % BmatRED = EIFEoper.INTforces.BmatRED(istrain,:);
        strainG = strainE(istrain) ;
        celas = EIFEoper.INTforces.MATPRO.celasglo(istrain,:) ;
        stressG = (celas*strainG) ;
        stressesREF{e}(:,g) = stressG ; 
        % ROTATION
        if  nstrain == 4
            ROTATION =  TRANSF_COORD{e}.ROTATION ;
            TRANSF_COORD{e}.ROTATION_STRESSES =    RotateStress2Dplanestrain(ROTATION(1,1),ROTATION(2,1)) ;
            stressG(1:3) =  TRANSF_COORD{e}.ROTATION_STRESSES*stressG(1:3) ;
        elseif nstrain == 6
            warning('Option not implemeted yet')
        else
            error('Option not implemented')
        end
        
        
        stressE(:,g) = stressG;
        stressAVERAGE  = stressAVERAGE + weig(g)*stressG ;
    end
    stressAVERAGEglo{e} = stressAVERAGE/sum(weig) ;
    stressN = sqrt(sum(stressE.^2,1)) ;
    [~,imax] = max(stressN) ;
    stressMAXglo{e} = stressE(:,imax) ;
    
    stressesROTATED{e} = stressE ; 
    
    if length(stressAVERAGEglo{e}) == 4
        % PLANE STRAIN
        % GID FORMAT IS DIFFERENT FROM OUR FORMAT
        %    MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy'}  ;
        stressAVERAGEglo{e} = stressAVERAGEglo{e}([1,2,4,3]) ;
        stressMAXglo{e} = stressMAXglo{e}([1,2,4,3]) ;
    elseif length(stressAVERAGEglo{e}) == 6
        stressAVERAGEglo{e} = stressAVERAGEglo{e}([1 2 3 6 4 5]) ;
        stressMAXglo{e} = stressMAXglo{e}([1 2 3 6 4 5]) ;
    else
        error('Option not implemented')
    end
    
    
    %
    %     stressGID_elem = zeros(nstrain,1) ;
    %     strainGID_elem = zeros(nstrain,1) ;
    %     weightELEM  = 0 ;
    %     for  g = 1:ngaus
    %         % Matrix of derivatives for Gauss point "g"
    %         BeXi = dershapef(:,:,g) ;
    %         % Jacobian Matrix
    %         Je = Xe*BeXi' ;
    %         detJe = det(Je) ;
    %         % Matrix of derivatives with respect to physical coordinates
    %         BeTILDE = inv(Je)'*BeXi ;
    %         % Matrix of symmetric gradient
    %         Be = QtransfB(BeTILDE,ndim) ;
    %         %
    %         strain = Be*dE;
    %         %  stress = celas*strain ;
    %         % Additional component
    %         if TP == 1
    %             strain3D = [strain(1) strain(2) 0 0 0 strain(3)]'; % Plane strain
    %             stress3D = celas3D*strain3D ;
    %             stressGID = stress3D([1 2 6 3]) ;  % For post-processing with GID
    %             strainGID = strain3D([1 2 6 3]) ;
    %         elseif TP == 2
    %             stress = celas*strain ;
    %             stress3D = [stress(1) stress(2) 0 0 0 stress(3)]'; % Plane strain
    %             strain3D = celas3Dinv*stress3D ;
    %             strainGID = strain3D([1 2 6 3]) ;
    %             stressGID = stress3D([1 2 6 3]) ; % For post-processing with GID
    %         elseif TP ==3
    %             stress = celas*strain ;
    %             stressGID = stress([1 2 3 6 4 5]) ; % For post-processing with GID
    %             strainGID = strain([1 2 3 6 4 5]) ; % For post-processing with GID
    %         end
    %         weight = detJe*weig(g);
    %         weightELEM = weightELEM + weight ;
    %         stressGID_elem = stressGID_elem + stressGID*weight ;
    %         strainGID_elem = strainGID_elem + strainGID*weight ;
    %
    %         indINI = (g-1)*nstrain+1 ;  indFIN = nstrain*g;
    %         stressGLO(indINI:indFIN,e) =  stressGID ;
    %         strainGLO(indINI:indFIN,e) =  strainGID ;
    %
end

stressAVERAGEglo = cell2mat(stressAVERAGEglo');
stressMAXglo = cell2mat(stressMAXglo');

stressDATA.MAXIM = stressMAXglo ;
stressDATA.AVERAGE = stressAVERAGEglo ;


USE_ROTATED_STRESSES_RECONST = 0; 
if USE_ROTATED_STRESSES_RECONST == 1
    stressesOUT = stressesROTATED ; 
else
     stressesOUT = stressesREF ;
end


%          stressGLO_elem(:,e) =  stressGID_elem/weightELEM ;
%         strainGLO_elem(:,e) =  strainGID_elem/weightELEM ;
%end
%DATA =DefaultField(DATA,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',0) ;
