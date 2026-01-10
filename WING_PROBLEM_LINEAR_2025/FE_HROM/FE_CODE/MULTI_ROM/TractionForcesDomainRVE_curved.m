function [FtractionE,Rast,rRB,fextDOMred,FORCES_2_PRINT_local ]= ...
    TractionForcesDomainRVE_curved(UPSILONforceTRAC,LOADS,...
    COOR,nnodeE,ndim,Q_tractions,...
    BETA_tractions,D_tractions,nmodes,ISLOCAL,DATAIN)

if nargin == 0
    load('tmp.mat')
end


% Vector of external forces due to traction loads
% -----------------------------------------------
FtractionE = zeros(nnodeE*max(ndim),1) ;
nfaces = min(length(UPSILONforceTRAC.COOR), length(LOADS));
Rast = zeros(nmodes.REACTION,1) ;

nRB = DATAIN.nRB ;
if nRB ==3
    ndimSP = 2;
else
    ndimSP = 3;
end


rRB =zeros(nRB,1) ;
fextDOMred =zeros(nmodes.DISPLACEMENT,1)  ;
FORCES_2_PRINT_local = cell(1,nfaces) ;
ndimGLO = ndim ;
ndim = max(ndim) ;

for iface =1:nfaces  % loop over all faces of the domain
    
    LOADF = LOADS{iface} ;
    
    if  ~isempty(LOADF)
        %  nmodes=  size(BETA_tractions{iface},1) ;
        % if nmodes >0
        % Rast = Rast + zeros(nmodes,1) ;
        %  fextDOMred = fextDOMred +  zeros(nmodes,1)  ;
        %end
        %   else
        F_face = cell(nnodeE,1) ;
        F_face(:) = {zeros(ndim,1) } ;
        COOR_DATA =  COOR{iface} ;  % Coordinates on which the input tractions data are defined
        COOR_FACE = UPSILONforceTRAC.COOR{iface}' ; % Coordinates of surface nodes
        
        %    ndimSP = 3; % Number of spatial dimensions
        nGAUSSbnd = size(UPSILONforceTRAC.OPER{1,iface},2)/ndimSP ; % Number of total Gauss points
        NORMALS = UPSILONforceTRAC.NORMALS{iface} ; % 3 x nelems. Unit normal vectors at each element of the face
        ngausELEM = nGAUSSbnd/size(NORMALS,2) ;  % number of Gauss points per element
        nELEMSS = nGAUSSbnd/ngausELEM ;  % number of elements
        
        
        if isempty(COOR_DATA) %|| (iscell(COOR_DATA) & isempty(COOR_DATA{1}))
            % It is assumed the input data is a uniform load
            
            % Therefore
             
              
              
            if ISLOCAL(iface) == 0
                
                 error('REvise this part. Introduce ROTATION MATRICES  !!!! ')
                 
                 
%                 % Rotation (future versions should perform this operation in a vectorized fashion)
%                 Tnodes = zeros(size(Tnodes_loc,1),nelem)  ; ;
%                 for ielem = 1:nelem
%                     Tnodes_loc =ROTATIONS_DOMAIN_transpose{ielem}*reshape(Tnodes_loc,ndimSP,[]) ;
%                     Tnodes(:,ielem)  = Tnodes_loc(:) ;
%                 end
%                 Tnodes = Tnodes(:) ; 
                
                Tuniform = LOADF(:) ;
                FORCES_2_PRINT_local{iface}  = repmat(Tuniform,nELEMSS,1) ;
                Tnodes = repmat(Tuniform,nGAUSSbnd,1) ;
                
            else
                % Forces aligned along the normal
                if ~isempty(find(LOADF(2:end)~=0))  
                    error('When ISLOCAL = 1, the nonzero component should be the first component')
                end
                
                Tnodes = zeros(ndimSP,nGAUSSbnd) ;
                for igaus = 1:ngausELEM
                    for idim = 1:ndimSP
                        Tnodes(idim,igaus:ngausELEM:end) = LOADF(1)*NORMALS(idim,:) ;
                    end
                end
                FORCES_2_PRINT_local{iface} =  Tnodes(:,1:ngausELEM:end) ;  % Force for each element
                Tnodes = Tnodes(:) ;
                
                
                
            end
            %      Tnodes = Tnodes(:) ;
        else
            
            
            %             Tnodes = zeros(ndimSP,nGAUSSbnd) ;
            %             for igaus = 1:ngausELEM
            %                 for idim = 1:ndimSP
            %                     Tnodes(idim,igaus:ngausELEM:end) = PRESION*NORMALS(idim,:) ;
            %                 end
            %             end
            %             FORCES_2_PRINT_local{iface} =  Tnodes(:,1:ngausELEM:end) ;  % Force for each element
            %             Tnodes = Tnodes(:) ;
            
            error('Program this part of the code !!!!')
            
        end
        
        for innodeE  = 1:nnodeE  % loop over interface  of domain (2 for slices)
            OPER = UPSILONforceTRAC.OPER{innodeE,iface} ; % Operator relating FE traction forces with beam forces
            F_face{innodeE} = OPER*Tnodes ;
        end
        
        Rast = Rast + Q_tractions{iface}*Tnodes ;
        fextDOMred = fextDOMred + D_tractions{iface}*Tnodes ;
        rRB = rRB + BETA_tractions{iface}*Tnodes ;
        
        
        if  any(ndimGLO(1)-ndimGLO)
            F_face_loc = cell(size(F_face)) ;
            for iloc = 1:length(F_face)
                nrows  = size(F_face{iloc}) ;
                F_face_loc{iloc} = zeros(ndim,1) ;
                F_face_loc{iloc}(1:nrows) = F_face{iloc} ;
                
            end
            FtractionE = FtractionE + cell2mat(F_face_loc) ;
            
        else
            FtractionE = FtractionE + cell2mat(F_face) ;
        end
        
        
        
    end
end
 