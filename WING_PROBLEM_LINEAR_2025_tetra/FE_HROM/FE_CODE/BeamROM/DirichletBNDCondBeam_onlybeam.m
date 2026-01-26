function [DOFr,DOFl,dR] = DirichletBNDCondBeam_onlybeam(DATAROM,MESH1D,DISP,DATAIN,ndim)

if nargin == 0
    load('tmp1.mat')
end

% JAHO, 25-Dec-2018 --- This option should set to zero always
DATAIN = DefaultField(DATAIN,'SET_TO_ZERO_NONRIGID_INTERFACE_MODES_AT_THE_BOUNDARIES',0) ; 


nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
ndimSP = size(MESH1D.COOR,2) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)



DISP = DefaultField(DISP,'LEFT_END',[]) ;
if ~isempty(DISP.LEFT_END)
    
    DISP.NODES{1} = DISP.LEFT_END  ;
    DISP.NODES{2} = DISP.RIGHT_END  ;
    MESH1D.NODES_POINTS{1} = MESH1D.LEFT_END_NODE ;
    MESH1D.NODES_POINTS{2} = MESH1D.RIGHT_END_NODE ;
    
    
end


if ndimSP == 2
    nrigid  =3 ;
else
    nrigid = 6;
end
% LEFT END
% ---------

% Known DOFs
dR = [] ;
DOFr = [] ;
  DATAROM = DefaultField(DATAROM,'nBOUNDARY_INTFMODES',[]) ;
  
%DATAIN = DefaultField(DATAIN,'ConsistentBoundaryConditions',1) ;   

for iend = 1:length(DISP.NODES)
    
    
    % This first block is only for rigid body modes 
    % ---------------------------------------------
    NODE{iend} = MESH1D.NODES_POINTS{iend} ;
    DOFS{iend} = small2large(NODE{iend},ndim) ;
    r = [] ;
    DISPLOC  =DISP.NODES{iend} ;
    for idim = 1:length(DISPLOC)
        if ~isempty(DISPLOC{idim})
            r = [r; idim] ;
            dR = [dR ; DISPLOC{idim}] ;
        end
    end
    DOFr = [DOFr; DOFS{iend}(r)];  ;
  
    % Not used 
     nbasisRIGID_BODY = nrigid ;
%     nBOUNDARY_INTFMODES = DATAROM.nBOUNDARY_INTFMODES ; % Number of interface boundary modes
%     if isempty(nBOUNDARY_INTFMODES)
%         nBOUNDARY_INTFMODES = 0 ;
%     end
    nBOUNDARY_INTFMODES = 0 ; 
    
    % Here we decide what to do with the remaining interface modes 
    % --------------------------------------------------------------
    ISSS =cellfun(@isempty,DISPLOC) ;
    if all(ISSS==0)   % All are prescribed displacements. We set accordingly the remaining DOFs to zero
        if  length(DOFS{iend}) >nbasisRIGID_BODY
            DOFr = [DOFr ; DOFS{iend}((nbasisRIGID_BODY+1):end)] ;
            dR = [dR; zeros(length(DOFS{iend}((nbasisRIGID_BODY+1):end)),1)] ; % Remaining DOFs
        end
    else
        % Free end (at least one DOFS)
        if   DATAIN.SET_TO_ZERO_NONRIGID_INTERFACE_MODES_AT_THE_BOUNDARIES  == 1   
            % We set to zero the amplitude of the non-rigid body modes at
            % the interfaces
            if  length(DOFS{iend}) >nbasisRIGID_BODY
                DOFr = [DOFr ;  DOFS{iend}((nbasisRIGID_BODY+1+nBOUNDARY_INTFMODES):end)] ;
                dR = [dR; zeros(length(DOFS{iend}((nbasisRIGID_BODY+1+nBOUNDARY_INTFMODES):end)),1)] ; % Remaining DOFs
            end
        end
    end
    
    
end






DOFl = setdiff(1:ndim*nnode,DOFr) ;



%
% % RIGHT END
% %----------------
%
%
% NODE_2 = MESH1D.RIGHT_END_NODE ;
% DOFS_2= small2large(NODE_2,ndim) ;
% % Known DOFs
% r = [] ;
% DISPLOC  =DISP.RIGHT_END ;
% for idim = 1:length(DISPLOC)
%     if ~isempty(DISPLOC{idim})
%         r = [r; idim] ;
%         dR = [dR ; DISPLOC{idim}] ;
%     end
% end
% DOFr = [DOFr; DOFS_2(r)] ;
%
%
%
% ISSS =cellfun(@isempty,DISPLOC) ;
%
% if all(ISSS==0)
%     if    length(DOFS_2) >6
%         DOFr = [DOFr ; DOFS_2((nbasisRIGID_BODY+1):end)] ;
%         dR = [dR; zeros(length(DOFS_2((nbasisRIGID_BODY+1):end)),1)] ;  % Remaining DOFs
%     end
% else
%     % Free end (at least one DOFS)
%     if  length(DOFS_2) >nbasisRIGID_BODY
%         DOFr = [DOFr ; DOFS_2((nbasisRIGID_BODY+1+nBOUNDARY_INTFMODES):end)] ;
%         dR = [dR; zeros(length(DOFS_2((nbasisRIGID_BODY+1+nBOUNDARY_INTFMODES):end)),1)] ; % Remaining DOFs
%     end
% end






