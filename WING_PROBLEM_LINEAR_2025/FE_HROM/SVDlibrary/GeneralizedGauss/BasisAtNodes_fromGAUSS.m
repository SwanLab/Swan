function [VAR_SMOOTH_FE,MSG,DATAIN] =  BasisAtNodes_fromGAUSS(DATA_GENGAUSS,NAME_INPUT_DATA,...
    BasisS_gauss,BdomRED_gauss,DATAIN,...
    DATA_REFMESH,HROMVAR,MSG,SingVal_F,VrightVal_F,wSTs)

if nargin == 0
    load('tmp1.mat')
end

% Checking incompatibilities
if DATAIN.TOL_LOC_InternalForces_IS_INTEGRATION == 1 || DATAIN.CUBATURE.IMPOSE_VOLUME_CONSTRAINT == 1
    error('Incompatibe options... ')
end
% FINDING NODAL EQUIVALENTS OF BasisS and BdomRED.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BasisS
% -------------------
ndim = size(DATA_REFMESH.COOR,2) ;
switch DATA_GENGAUSS.APPROX_FUN__DERI.METHOD
    case {'FE_INTERPOLATION'}
        DATAINloc.ncomponent = HROMVAR.nstrain ;
        DATAINloc.RETURN_GAUSS_APPROX = 1;
        [VARNODES,VARGAUSS]= GaussToNodesSmooth([BasisS_gauss,BdomRED_gauss],DATA_REFMESH.Nst(1:ndim:end,1:ndim:end),DATAINloc) ;
        
        IND = 1:size(BasisS_gauss,2) ;
        BasisS_nodes = VARNODES(:,IND) ;
        BasisS_gauss_approx = VARGAUSS(:,IND) ;
        MSG =  dispMSG(MSG, '------------------------------------------')  ;
        MSG =  dispMSG(MSG, 'FE-based interpolation method')  ;
        MSG  = dispMSG(MSG, '------------------------------------------')  ;
        MSG  = dispMSG(MSG, '------------------------------------------')  ;
        ERROR_intBasisS = BasisS_gauss_approx-BasisS_gauss ;
        nERROR = norm(ERROR_intBasisS,'fro')/norm(BasisS_gauss,'fro') ;
        MSG  = dispMSG(MSG, ['BasisS --> FROB. ERROR transfer. Gauss-to-nodes (%) =',num2str(nERROR*100)])  ;
        
        
        IND = size(BasisS_gauss,2)+1;  ;
        BdomRED_nodes = VARNODES(:,IND:end) ;
        BdomRED_gauss_approx = VARGAUSS(:,IND:end) ;
        
        ERROR_intB = BdomRED_gauss_approx-BdomRED_gauss ;
        nERROR = norm(ERROR_intB,'fro')/norm(BdomRED_gauss,'fro') ;
        MSG  = dispMSG(MSG, ['BdomRED --> FROB. ERROR transfer. Gauss-to-nodes (%) =',num2str(nERROR*100)])  ; ,
        MSG  = dispMSG(MSG, '------------------------------------------')  ;
        
        
        
        DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_ERROR_IN_SMOOTHING_Fint',1) ;
        
        if DATA_GENGAUSS.PLOT_ERROR_IN_SMOOTHING_Fint == 1
            nstrain =  HROMVAR.nstrain ;
            ngaus = size(BasisS_gauss,1)/nstrain ;
            SNAPforceS = BasisRED_BasisStress(BdomRED_gauss,BasisS_gauss,nstrain,ngaus) ;  % Original Internal Force MAtrix
            SNAPforceS_approx = BasisRED_BasisStress(BdomRED_gauss_approx,BasisS_gauss_approx,nstrain,ngaus) ;
            
            ERROR_intF = SNAPforceS_approx-SNAPforceS ;
            nERROR = norm(ERROR_intF,'fro')/norm(SNAPforceS,'fro') ;
            MSG  = dispMSG(MSG, ['Internal Forces --> FROB. ERROR transfer. Gauss-to-nodes (%) =',num2str(nERROR*100)])  ;
            MSG  = dispMSG(MSG, '------------------------------------------')  ;
            DATAIN.LABEL_NAME_VARIABLE = 'ErrSmoothFint(t.p.)' ;
            DATAIN.LABEL_NAME_FILE = 'ErrSmF' ;
            % ERROR_intF = ERROR_intF*VAR_SMOOTH_FE.invSVsingular_F   ;% DOES NOT WORK !  To get the error in approximating the modes
            
            % Let us display ERROR_intF in a normalized way . For each snapshot, we calculate its average value
            
            for imode = 1:size(ERROR_intF,2)
                % Reference value
             refVAL = max(abs(SNAPforceS(:,imode))) ; 
              %  refVAL = sqrt(sum(refVAL.^2))/sum(wSTs) ;
              %  ERROR_intF(:,imode) = ERROR_intF(:,imode)./SNAPforceS(:,imode)*100 ;
                
                  ERROR_intF(:,imode) = ERROR_intF(:,imode)./refVAL*100 ;
            end
            
            DATAIN =   GidPostProcessModesFINT_loc(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
                DATA_REFMESH.TypeElement,ERROR_intF,DATA_REFMESH.posgp,NAME_INPUT_DATA,DATAIN);
            
            for iii = 1:length(DATAIN.MSGPRINT)
                MSG{end+1} = DATAIN.MSGPRINT{iii} ;
            end
            
        end
        
        
        
        % Variables for performing the spatial search
        % -----------------------------------------------
        % Reduced B-matrix at nodes
        VAR_SMOOTH_FE.BdomRED_nodes = BdomRED_nodes;
        % Stress basis matrix at nodes
        VAR_SMOOTH_FE.BasisS_nodes = BasisS_nodes;
        
    case 'FE_ONLY'
        
        
        MSG  = dispMSG(MSG, '------------------------------------------')  ;
        MSG  = dispMSG(MSG, ['MEthod with interpolation using Gauss points as interpolation nodes']  );
        MSG  = dispMSG(MSG, '------------------------------------------')  ;
        
        % Reduced B-matrix at nodes
        VAR_SMOOTH_FE.BdomRED_gauss = BdomRED_gauss;
        % Stress basis matrix at nodes
        VAR_SMOOTH_FE.BasisS_gauss = BasisS_gauss;
        
        
         VAR_SMOOTH_FE.BdomRED_nodes = [];
        % Stress basis matrix at nodes
        VAR_SMOOTH_FE.BasisS_nodes = [];
        
      
        
end


VAR_SMOOTH_FE.ngausE = size(DATA_REFMESH.posgp,2); 
% Inverse of the singular values of BasisF times its right-singular vectors
% --------------------------------------------------------------------------
VAR_SMOOTH_FE.invSVsingular_F = bsxfun(@times,VrightVal_F',1./SingVal_F)' ;




% Matrix of connectivities
% ---------------------------
VAR_SMOOTH_FE.COOR = DATA_REFMESH.COOR ;
% Matrix of coordinates
% ---------------------------
VAR_SMOOTH_FE.CN = DATA_REFMESH.CN ;
% Type of element
% --------------------------
VAR_SMOOTH_FE.TypeElement = DATA_REFMESH.TypeElement ;
% Dela. triangulation
% --------------------------
VAR_SMOOTH_FE.DELTRIANG = delaunayTriangulation(DATA_REFMESH.COOR);
% ORDER POLYNOMIALS
% ----------------------------Ç
nnodeE = size(DATA_REFMESH.CN,2) ;
switch VAR_SMOOTH_FE.TypeElement
    case 'Quadrilateral'
        IND_POLYG_ELEMENT = [1 2 3 4 1] ;
        if nnodeE == 9
            ORDER_POLYNOMIALS  =[2 2] ;
            
        elseif nnodeE == 4
            ORDER_POLYNOMIALS  =[1 1] ;
        else
            error('element not implemented')
        end
    case 'Hexahedra'
        IND_POLYG_ELEMENT = [1 2 3 4 5 6 7 8 1] ;
        if nnodeE == 27
            ORDER_POLYNOMIALS  =[2 2 2] ;
        elseif nnodeE == 8
            ORDER_POLYNOMIALS  =[1 1 1] ;
        else
            error('element not implemented')
        end
    otherwise
        error('element not implemented')
end
VAR_SMOOTH_FE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
VAR_SMOOTH_FE.IND_POLYG_ELEMENT = IND_POLYG_ELEMENT;   % Local numbering of corner nodes (polygon)
VAR_SMOOTH_FE.nstrain = HROMVAR.nstrain  ;


