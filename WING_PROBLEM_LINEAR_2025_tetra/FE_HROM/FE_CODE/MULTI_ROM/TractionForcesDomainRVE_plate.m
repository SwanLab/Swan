function [FtractionE,Rast,rRB,fextDOMred,FORCES_2_PRINT_local ]= TractionForcesDomainRVE_plate(UPSILONforceTRAC,LOADS,...
    COOR,nnodeE,ndim,Q_tractions,...
    BETA_tractions,D_tractions,nmodes,ISLOCAL)




% Vector of external forces due to traction loads
% -----------------------------------------------
FtractionE = zeros(nnodeE*ndim,1) ;
nfaces = min(length(UPSILONforceTRAC.COOR), length(LOADS));
Rast = zeros(nmodes.REACTION,1) ;
nmodesRB = 6 ;
rRB =zeros(nmodesRB,1) ;
fextDOMred =zeros(nmodes.DISPLACEMENT,1)  ;
FORCES_2_PRINT_local = cell(1,nfaces) ;

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
        
        ndimSP = 3; % Number of spatial dimensions
        nGAUSSbnd = size(UPSILONforceTRAC.OPER{1,iface},2)/ndimSP ; % Number of total Gauss points
        NORMALS = UPSILONforceTRAC.NORMALS{iface} ; % 3 x nelems. Unit normal vectors at each element of the face
        ngausELEM = nGAUSSbnd/size(NORMALS,2) ;  % number of Gauss points per element
        nELEMSS = nGAUSSbnd/ngausELEM ;  % number of elements
        
        
        if isempty(COOR_DATA) %|| (iscell(COOR_DATA) & isempty(COOR_DATA{1}))
            % It is assumed the input data is a uniform load
            
            % Therefore
            
            if ISLOCAL(iface) == 0
                Tuniform = LOADF(:) ;
                FORCES_2_PRINT_local{iface}  = repmat(Tuniform,nELEMSS,1) ;
                Tnodes = repmat(Tuniform,nGAUSSbnd,1) ;
                
            else
                % Forces aligned along the normal
                if LOADF(1) == 0
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
        
        %for innodeE  = 1:nnodeE  % loop over interface  of domain (2 for slices)
%             OPER = UPSILONforceTRAC.OPER{innodeE,iface} ; % Operator relating FE traction forces with beam forces
%             F_face{innodeE} = OPER*Tnodes ;
             OPER = UPSILONforceTRAC.OPER{iface} ; % Operator relating FE traction forces with beam forces
            F_face  = OPER*Tnodes ;
        %end
        
        Rast = Rast + Q_tractions{iface}*Tnodes ;
        fextDOMred = fextDOMred + D_tractions{iface}*Tnodes ;
        rRB = rRB + BETA_tractions{iface}*Tnodes ;
        
        FtractionE = FtractionE +  (F_face) ;
        
    end
end
