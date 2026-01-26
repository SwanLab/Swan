function [TRANS_COOR,normD] = PolarDecompEIFEperm(Xe,EIFEoper,DATA)
% This function determines the rotation of the fictitious interface modes given the coordinates
% of the EIF element Xe, as well as the properties of the parent element
% contained in the data structure % EIFEoper
% JAHO,2-MArch-2023
% ----------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
elseif  nargin == 2
    DATA = [] ; 
    DATA.ToleranceDeformationalPartTransformationParentDomain = 1e-4 ; 
end
[ndim,nnodeE] = size(Xe) ;
 
% COMPUTATION OF THE DILATION FACTOR (ACTUALLY THE DETERMINANT OF THE
% JACOBIAN)
% AS WELL AS THE ROTATION MATRIX SO THAT THE DILATATIONAL PART OF
% THE DEFORMATION IS ZERO (AND IF CANNOT BE RENDER ZERO, THEN THE PROGRAM
% WILL COME TO A HALT, BECAUSE THE CURRENT VERSION (3-March-2023) cannot
% handle accurately distorted domains )
TRANS_COOR = [] ;
METHOD = 2;  % METHOD = 1, METHOD = 2, both methods yield identical results
if METHOD == 1
    
    % THIS METHOD employs the JACOBIAN OF THE TRANSFORMATION
    % METHOD 2 DOESN'T, AND RATHER RESORTS TO MINIMIZATION OF THE
    % DILATATIONAL PART
    % -----------------------------------------------------
    
    % Jacobian matrix at point g  (transformation of coordinates )
    % Matrix of derivatives for Gauss point "g" (polynomial shape functions)
    %  BeXi = dershapef(:,:,g) ;
    g = 1;
    BeXi_transf = zeros(ndim,nnodeE) ;
    for idim = 1:ndim
        BeXi_transf(idim,:) =  EIFEoper.INTforces.Bgrad_transf{idim}(g,:) ;
    end
    % Jacobian Matrix
    %     c = sum(Xe,1)/3  ;
    %     for idim = 1:3
    %         Xe(:,idim) = Xe(:,idim)-c(idim);
    %     end
    
    Je = Xe*BeXi_transf' ;
    % JAcobian
    detJe = det(Je) ;
    
    [UU,SS,VV] = SVDT(Je) ;
    % Decomposition  Je = Q*D
    Q = UU*VV' ;  % Rotational part
    if abs(norm( (SS)/SS(1)) - sqrt(size(Je,2))) > 1e-10
        % disp('Dilatational part different from zero')
        % Vall_rot = [] ;
    else
        
        
        %         Vall_rot = zeros(size(EIFEoper.MODES.Vall)) ;
        %         for imodes = 1:size(EIFEoper.MODES.Vall,2)
        %             LOCV = (Je'\reshape(EIFEoper.MODES.Vall(:,imodes),ndim,[])) ;
        %             %   LOCV = (Q*reshape(EIFEoper.MODES.Vall(:,imodes),ndim,[]))/detJe ;
        %             Vall_rot(:,imodes) = LOCV(:) ;
        %         end
        
        
        % TRANSLATION
        TRANSLATION_ce = sum(Xe,2)/size(Xe,2) ;
        
        TRANS_COOR.ROTATION = Q;
        TRANS_COOR.TRANSLATION = TRANSLATION_ce;
        TRANS_COOR.SCALEFACTOR = (detJe)^(1/ndim);
        TRANS_COOR.detJe = detJe ;
        
        %    Xe = TRANSLATION_ce + TRANS_COOR.SCALEFACTOR*Q*EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord.COOR';
        
    end
    
    
    
    
else
    
    % THIS METHOD DOES NOT REQUIRE TO COMPUTE THE JACOBIAN OF THE
    % TRANSFORMATION
    
    
    % X_e = TRANSLATION_ce + DILAT_FACTOR*Q*(Xref) ;
    % The preceding expression assumes that Xref are the coordinates of the
    % parent domain with respect to its centroid.
    % However, such centroid is calculated using the underlying FE mesh. In
    % the case of the physical domain, we cannot calculate such centroid
    % because there are no underlying FE mesh.
    % For instance, if there is neither rotation nor dilation, then
    % X_e  -Ce  =    Xref - C_ref ---> X_e = C_e + (Xref-C_ref)
    % If we introduce rotations
    % X_e = C_e + DILAT_FACTOR*Q*(Xref-C_ref) ;
    
    
    
    Xref=  EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord.COOR'; % Coordinates parent domain
    CNref = EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord.CN;  % Connectivities parent domain 
    
    PLOT_xref = 0; 
    if PLOT_xref == 1
        figure(4)
        hold on
        for iplott = 1:size(Xref,2)
            plot(Xref(1,iplott),Xref(2,iplott)) ; 
            text(Xref(1,iplott),Xref(2,iplott),num2str(iplott))
        end
        
    end
    
     DATA = DefaultField(DATA,'MESH_GIVEN_BY_THE_USER',[]) ; 
    Xref= Xref(:,CNref) ;
    if (size(CNref,2) == 27 || size(CNref,2) == 9) && isempty(DATA.MESH_GIVEN_BY_THE_USER)
      
        Xref = Xref(:,1:end-1) ;
   
    end
    
    
    
    C_ref = sum(Xref,2)/size(Xref,2) ; % Pseudo-centroid parent domain 
    
    XrefREL = bsxfun(@minus,Xref,C_ref) ; % Relative coordinates parent domain 
    
    % Compute maximum length parent domain 
    % JAHO, 17-APRIL-2024
    % It is assumed here that the connectivities are the same for both the
    % reference domain and the pyhysical domain 
    % --------------------------------
    METHOD_COMPUTE_maximum_length = 0;  % 0  Before 17-Apr-2024
    dLref = DistanceMaximPointsEIFEM(Xref,METHOD_COMPUTE_maximum_length) ; 
    % -------------------------------------
    % Translation. First compuyte the pseudo-centroid of the physical
    % domain 
    C_e = sum(Xe,2)/size(Xe,2) ;
    XeREL = bsxfun(@minus,Xe,C_e) ; % Relative coordinates
    TRANSLATION_ce = C_e ;  % This is the translation 
    % --------------------------------------
    % Scale factor
%     XAUG = [XeREL,XeREL(:,1)] ;
%     dL  =diff(XAUG,1,2) ;
%     dLE = max(sqrt(sum(dL.^2,1))) ;
     dLE = DistanceMaximPointsEIFEM(XeREL,METHOD_COMPUTE_maximum_length) ; 
    factorSCALE = dLref/dLE ;
    %--------------------------------------
    XeRELs = XeREL*factorSCALE;
    % ---------------------------------------
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/05_linearTRIANGLE.mlx
    % --------------------------------------
    [UU,SS,VV]  =SVDT(XrefREL*XeRELs') ;
    
     % Ensure that the determinant of Q is 1 (not -1)
    Sdet = eye(ndim) ;
    Sdet(ndim,ndim) = det(UU*VV') ;
    Q = (UU*Sdet*VV')' ;
    % Deformational part
    D = XeRELs -Q*XrefREL;
    TOL = DATA.ToleranceDeformationalPartTransformationParentDomain; % = DefaultField(DATA,'',1e-5) ;
    normD = norm(D,'fro')/length(D); 
  %  disp(['Norm def. part = ',num2str(normD)])
   %   normDmax = min(normDmax,normD) ; 
    
    PLOTTR = 0 ;
    %close all
    PermutSearchParentPlot(XeRELs,XrefREL,PLOTTR,Q'*XeRELs) ; 
    
   

    % disp(['Norm deformational part = ',num2str(norm(D,'fro')/length(D))])
    if normD >= TOL
        %  error('Deformational part must be negligible')
        %Vall_rot = [] ;
        %detJe = [] ;
      
    else
        
        %         Vall_rot = zeros(size(EIFEoper.MODES.Vall)) ;
        %         for imodes = 1:size(EIFEoper.MODES.Vall,2)
        %             LOCV = factorSCALE*(Q*reshape(EIFEoper.MODES.Vall(:,imodes),ndim,[])) ;
        %             %   LOCV = (Q*reshape(EIFEoper.MODES.Vall(:,imodes),ndim,[]))/detJe ;
        %             Vall_rot(:,imodes) = LOCV(:) ;
        %         end
        detJe = 1/(factorSCALE^ndim) ;
        TRANS_COOR.ROTATION = Q;
        TRANS_COOR.TRANSLATION = TRANSLATION_ce;
        TRANS_COOR.CenterRotationREFERENCE = C_ref;
        TRANS_COOR.SCALEFACTOR = 1/factorSCALE;
        TRANS_COOR.detJe = detJe ;
        %  Xe = TRANSLATION_ce + 1/factorSCALE*Q*Xref;
        
    end
    
    
    
    
    
    
    
    
    
end



%     else
%         if abs(detJe-detJe_ref)/detJe_ref >= 1e-10
%             error('Only constant Jacobian transformations allowed')
%         end