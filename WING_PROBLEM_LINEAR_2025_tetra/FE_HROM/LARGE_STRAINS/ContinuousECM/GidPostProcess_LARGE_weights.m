function [NameFile_res ]= GidPostProcess_LARGE_weights(COOR,CN,TypeElement, ...
    posgp,NameFile_res,MaterialType,IND_ELEM_MESHES,...
    MESH,DATA,WEIGHTS);
% Post-processing of results using GID
%dbstop('5')
if nargin==0
    load('tmp1.mat')
end


% Writing results file
% -----------------------------------------------------
% Nodal vectors
% -----------------------------------------------------
NODAL_VECTOR = [] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------
% STRESS DATA. 1D mesh
% -----------------------



nelem_TOT =  0 ;
ielemINI = 1; 

for imesh = 1:(length(TypeElement)-1)
    MESH(imesh).TypeElement = TypeElement{imesh} ;
    MESH(imesh).posgp = posgp{imesh} ;
    MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
    nelemE(imesh) = size(CN{imesh},2) ;
    nelem_TOT = nelem_TOT +  size(CN{imesh},1) ;
    
    %if ~isempty(STRESSDATA.MAX_VONMISES)
    % scalar defined on GAUSS POINTS, mesh 1
    
    MESH(imesh).GAUSS_SCALAR(1).NAME = 'WEIGHTS' ;
    MESH(imesh).GAUSS_SCALAR(1).VAR =WEIGHTS{imesh} ;
    ielemFIN = ielemINI +size(CN{imesh},1)-1  ; 
    MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = ielemINI:ielemFIN  ; % IND_ELEM_MESHES{imesh} ;
    ielemINI = ielemFIN+1; 
end

%end

%
% %%%%%%%%%%%%%%%%%%%%%%%
% % 2D represetnations
% DATA = DefaultField(DATA,'GeneralizedForces',[]) ;
%
% if ~isempty(DATA.GeneralizedForces)
%     fff = fieldnames(DATA.GeneralizedForces) ;
%     for iresults =1:length(fff)
%         INFO = DATA.GeneralizedForces.(fff{iresults}) ;
%         MESH(imesh).GAUSS_VECTOR(iresults).NAME = INFO.NAME ;
%         MESH(imesh).GAUSS_VECTOR(iresults).COMP = INFO.LEGEND ;
%         MESH(imesh).GAUSS_VECTOR(iresults).VAR = INFO.VALUE ;
%         MESH(imesh).GAUSS_VECTOR(iresults).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%     end
% end



%
%
% RESIDUAL = DefaultField(RESIDUAL,'MAX_NORM',[]) ;
%
% if ~isempty(RESIDUAL.MAX_NORM)
%     if DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES > 1
%         APPLEG = ['(x',num2str(DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES),')'] ;
%     end
%     MESH(imesh).GAUSS_SCALAR(2).NAME = ['RESIDF_Max_Norm',APPLEG] ;
%     MESH(imesh).GAUSS_SCALAR(2).VAR = RESIDUAL.MAX_NORM ;
%     MESH(imesh).GAUSS_SCALAR(2).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%     MESH(imesh).GAUSS_SCALAR(3).NAME = ['RESIDF_Norm_Avg',APPLEG] ;
%     MESH(imesh).GAUSS_SCALAR(3).VAR = RESIDUAL.NORM_AVERAGE ;
%     MESH(imesh).GAUSS_SCALAR(3).ELEMENTS = IND_ELEM_MESHES{imesh} ;
% end

% ------------------------------------------------------
%%%%%%%% VARIABLES DEFINED ON 3D GAUSS POINT MESHES
% % -----------------------------------------------------
% STRESSDATA  = DefaultField(STRESSDATA,'VONMISES',[]) ;
%
% for imeshLOC  = 1:length(STRESSDATA.VONMISES)
%     imesh = imeshLOC + 1;   % The first mesh is 1D
%
%     if ~isempty(CN{imesh})
%         % ----------
%         % 3D mesh
%         % ----------
%         MESH(imesh).TypeElement = TypeElement{imesh} ;
%         MESH(imesh).posgp = posgp{imesh} ;
%         MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
%         nelemE(imesh) = size(CN{imesh},2) ;
%
%         iRESULTS = imesh -1 ;
%
%         if ~isempty(STRESSDATA.VONMISES)
%             MESH(imesh).GAUSS_SCALAR(1).NAME = 'Von_Mises_3D' ;
%             MESH(imesh).GAUSS_SCALAR(1).VAR = STRESSDATA.VONMISES{iRESULTS}' ;
%             MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%         end
%         %
%
%         if ~isempty(stressGLO)
%             MESH(imesh).GAUSS_MATRIX(1).NAME = 'STRESS' ;
%             MESH(imesh).GAUSS_MATRIX(1).VAR = stressGLO{iRESULTS} ;
%             MESH(imesh).GAUSS_MATRIX(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%          %   if size(stressGLO{iRESULTS},1) == 6
%                 MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
%
%           %  elseif  size(stressGLO{iRESULTS},1) == 4
%            %     MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy'}  ;
%            % else
%             %    error('Option not implemented')
%            % end
%
%         end
%         %  nelemE(2) = size(CN{2},2) ;
%
% %         DATAIN = DefaultField(DATAIN,'POST_PROCESS_TRACTION_FORCES',0) ; % = 1;  % New option, 30-July-2019
% %
% %         if DATAIN.POST_PROCESS_TRACTION_FORCES == 1
% %              [Tnod, ElemLatSurf]= TractionForcesPrint(MESH,FORCES_2_PRINT,DATAIN) ;
% %
% %                % External traction forces
% % %         % ----------------- %
% %         MESH(imesh).GAUSS_VECTOR(1).NAME = 'TRACTION' ;
% %          MESH(imesh).GAUSS_VECTOR(1).VAR   = Tnod ;
% %          MESH(imesh).GAUSS_VECTOR(1).COMP = {'t-x','t-y','t-z'} ;
% %          MESH(imesh).GAUSS_VECTOR(1).ELEMENTS = ElemLatSurf' + nelem_TOT  ;
% %
% %         end
%
% %
% %
% %
%
%
%     else
%         % Lateral surface mesh
%         % ----------------------
%         % ----------
%
%         imesh = 2;
%         MESH(imesh).TypeElement = TypeElement{3} ;
%         MESH(imesh).posgp = posgp{3} ;
%         MESH(imesh).NAMEMESH = DATA.NAMEMESH{3};
%         nelemE(imesh) = size(CN{3},2) ;
%
%
%         [Tnod, ElemLatSurf]= TractionForcesPrint(MESH,FORCES_2_PRINT,DATAIN) ;
%
%         % External traction forces
%         % ----------------- %
%         MESH(imesh).GAUSS_VECTOR(1).NAME = 'TRACTION' ;
%         MESH(imesh).GAUSS_VECTOR(1).VAR   = Tnod ;
%         MESH(imesh).GAUSS_VECTOR(1).COMP = {'t-x','t-y','t-z'} ;
%         MESH(imesh).GAUSS_VECTOR(1).ELEMENTS = ElemLatSurf' + nelem_TOT  ;
%
%
%     end
%
% end

ndim = size(COOR,2) ;
DATALOC = [] ;
DATAIN.TIMELOC = 1; % = DefaultField(DATAIN,'TIMELOC',1) ;
DATALOC.TIMELOC = DATAIN.TIMELOC;
NODAL_SCALAR = [] ;
GidResults2DFE_multi(NameFile_res,ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATALOC);
