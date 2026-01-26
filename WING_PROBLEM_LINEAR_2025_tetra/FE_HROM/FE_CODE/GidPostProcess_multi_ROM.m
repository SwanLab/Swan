function [NameFile_res ]= GidPostProcess_multi_ROM(COOR,CN,TypeElement,d,strainGLO, stressGLO, ...
    React,posgp,NameFile_res,MaterialType,NODES,RESIDUAL,DATAIN,DATA,STRESSDATA,IND_ELEM_MESHES,...
    FORCES_2_PRINT,DATA_REFMESH,DATARUN);
% Post-processing of results using GID
%dbstop('5')
if nargin==0
    load('tmp2.mat')
end


% Writing results file
% -----------------------------------------------------
% Nodal vectors
% -----------------------------------------------------
NODAL_VECTOR(1).NAME = 'DISPLACEMENTS' ;
NODAL_VECTOR(1).COMP = {'X-DISP','Y-DISP','Z-DISP'};
NODAL_VECTOR(1).VAR = d ;
NODAL_VECTOR(1).NODES = 1:size(COOR,1) ;

DATA = DefaultField(DATA,'REACTIONS3D',[]) ; 

if ~isempty(DATA.REACTIONS3D)
    NODAL_VECTOR(2).NAME = 'REACTIONS-SE' ;
    NODAL_VECTOR(2).COMP = {'X-REAC','Y-REAC','Z-REAC'};
    NODAL_VECTOR(2).VAR = DATA.REACTIONS3D ;
    NODAL_VECTOR(2).NODES = 1:size(COOR,1) ;
end



%
%
% if ~isempty(RESIDUAL.ALL)
%     if DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES > 1
%         APPLEG = ['(x',num2str(DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES),')'] ;
%     end
%
%
%     NODAL_VECTOR(2).NAME = ['RES.FORCES',APPLEG] ;
%     NODAL_VECTOR(2).COMP = {'X-RES','Y-RES','Z-RES'};
%     NODAL_VECTOR(2).VAR = RESIDUAL.ALL ;
%     NODAL_VECTOR(2).NODES = NODES{2}' ;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar vectors: forces and moments
% -------------------------------------
DATA = DefaultField(DATA,'Forces1D',[]) ;
NODAL_SCALAR = [] ; iacum = 1 ;
if ~isempty(DATA.Forces1D)
    DATAIN  = DefaultField(DATAIN,'FACTOR_MULTIPLY_MOMENTS_DIAGRAM',1) ;
    if DATAIN.FACTOR_MULTIPLY_MOMENTS_DIAGRAM > 1
        APPLEG = ['(x',num2str(DATAIN.FACTOR_MULTIPLY_MOMENTS_DIAGRAM),')'] ;
    else
        APPLEG ='';
    end
    fff = fieldnames(DATA.Forces1D) ;
    for ifield = 1:length(fff)
        FORCE = DATA.Forces1D.(fff{ifield}) ;
        if ~isempty(FORCE)
            NODAL_SCALAR(iacum).NAME = [fff{ifield},APPLEG ];
            NODAL_SCALAR(iacum).VAR = FORCE(:) ;
            NODAL_SCALAR(iacum).NODES = NODES{1} ;
            iacum = iacum + 1;
        end
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------
% STRESS DATA. 1D mesh
% -----------------------

STRESSDATA  = DefaultField(STRESSDATA,'MAX_VONMISES',[]) ;

nelem_TOT =  0 ;

imesh = 1 ;
MESH(imesh).TypeElement = TypeElement{imesh} ;
MESH(imesh).posgp = posgp{imesh} ;
MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
nelemE(imesh) = size(CN{imesh},2) ;
nelem_TOT = nelem_TOT +  size(CN{imesh},1) ;

if ~isempty(STRESSDATA.MAX_VONMISES)
    % scalar defined on GAUSS POINTS, mesh 1
    DATAIN = DefaultField(DATAIN,'POST_PROCESS_VON_MISES_PER_MATERIAL',0);  % 19-Feb-2020
    if DATAIN.POST_PROCESS_VON_MISES_PER_MATERIAL > 0
        EXTRL = [' MT=',num2str(DATAIN.POST_PROCESS_VON_MISES_PER_MATERIAL)] ; 
    else
        EXTRL = '' ; 
    end
    
    MESH(imesh).GAUSS_SCALAR(1).NAME = ['MaxVonMises',EXTRL] ;
    MESH(imesh).GAUSS_SCALAR(1).VAR = STRESSDATA.MAX_VONMISES ;
    MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%
% 2D represetnations
DATA = DefaultField(DATA,'GeneralizedForces',[]) ;

if ~isempty(DATA.GeneralizedForces)
    fff = fieldnames(DATA.GeneralizedForces) ;
    for iresults =1:length(fff)
        INFO = DATA.GeneralizedForces.(fff{iresults}) ;
        MESH(imesh).GAUSS_VECTOR(iresults).NAME = INFO.NAME ;
        MESH(imesh).GAUSS_VECTOR(iresults).COMP = INFO.LEGEND ;
        MESH(imesh).GAUSS_VECTOR(iresults).VAR = INFO.VALUE ;
        MESH(imesh).GAUSS_VECTOR(iresults).ELEMENTS = IND_ELEM_MESHES{imesh} ;
    end
end





RESIDUAL = DefaultField(RESIDUAL,'MAX_NORM',[]) ;

if ~isempty(RESIDUAL.MAX_NORM)
    if DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES > 1
        APPLEG = ['(x',num2str(DATAIN.FACTOR_MULTIPLY_RESIDUAL_FORCES),')'] ;
    end
    MESH(imesh).GAUSS_SCALAR(2).NAME = ['RESIDF_Max_Norm',APPLEG] ;
    MESH(imesh).GAUSS_SCALAR(2).VAR = RESIDUAL.MAX_NORM ;
    MESH(imesh).GAUSS_SCALAR(2).ELEMENTS = IND_ELEM_MESHES{imesh} ;
    MESH(imesh).GAUSS_SCALAR(3).NAME = ['RESIDF_Norm_Avg',APPLEG] ;
    MESH(imesh).GAUSS_SCALAR(3).VAR = RESIDUAL.NORM_AVERAGE ;
    MESH(imesh).GAUSS_SCALAR(3).ELEMENTS = IND_ELEM_MESHES{imesh} ;
end

% ------------------------------------------------------
%%%%%%%% VARIABLES DEFINED ON 3D GAUSS POINT MESHES
% -----------------------------------------------------
STRESSDATA  = DefaultField(STRESSDATA,'VONMISES',[]) ;

for imeshLOC  = 1:length(STRESSDATA.VONMISES)
    imesh = imeshLOC + 1;   % The first mesh is 1D
    
    if ~isempty(CN{imesh})
        % ----------
        % 3D mesh
        % ----------
        MESH(imesh).TypeElement = TypeElement{imesh} ;
        MESH(imesh).posgp = posgp{imesh} ;
        MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
        nelemE(imesh) = size(CN{imesh},2) ;
        
        iRESULTS = imesh -1 ;
        
        if ~isempty(STRESSDATA.VONMISES)
            MESH(imesh).GAUSS_SCALAR(1).NAME = 'Von_Mises_3D' ;
            MESH(imesh).GAUSS_SCALAR(1).VAR = STRESSDATA.VONMISES{iRESULTS}' ;
            MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
        end
        %
        
        if ~isempty(stressGLO)
            MESH(imesh).GAUSS_MATRIX(1).NAME = 'STRESS' ;
            MESH(imesh).GAUSS_MATRIX(1).VAR = stressGLO{iRESULTS} ;
            MESH(imesh).GAUSS_MATRIX(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
         %   if size(stressGLO{iRESULTS},1) == 6
                MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
                
          %  elseif  size(stressGLO{iRESULTS},1) == 4
           %     MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy'}  ;
           % else
            %    error('Option not implemented')
           % end
            
        end
        %  nelemE(2) = size(CN{2},2) ;
        
%         DATAIN = DefaultField(DATAIN,'POST_PROCESS_TRACTION_FORCES',0) ; % = 1;  % New option, 30-July-2019 
% 
%         if DATAIN.POST_PROCESS_TRACTION_FORCES == 1 
%              [Tnod, ElemLatSurf]= TractionForcesPrint(DATA_REFMESH,FORCES_2_PRINT,DATAIN) ;
%              
%                % External traction forces
% %         % ----------------- %
%         MESH(imesh).GAUSS_VECTOR(1).NAME = 'TRACTION' ;
%          MESH(imesh).GAUSS_VECTOR(1).VAR   = Tnod ;
%          MESH(imesh).GAUSS_VECTOR(1).COMP = {'t-x','t-y','t-z'} ;
%          MESH(imesh).GAUSS_VECTOR(1).ELEMENTS = ElemLatSurf' + nelem_TOT  ;
%              
%         end
        
%         
%         
%       
        
        
    else
        % Lateral surface mesh
        % ----------------------
        % ----------
        
        imesh = 2;
        MESH(imesh).TypeElement = TypeElement{3} ;
        MESH(imesh).posgp = posgp{3} ;
        MESH(imesh).NAMEMESH = DATA.NAMEMESH{3};
        nelemE(imesh) = size(CN{3},2) ;
        
        
        [Tnod, ElemLatSurf]= TractionForcesPrint(DATA_REFMESH,FORCES_2_PRINT,DATAIN) ;
        
        % External traction forces
        % ----------------- %
        MESH(imesh).GAUSS_VECTOR(1).NAME = 'TRACTION' ;
        MESH(imesh).GAUSS_VECTOR(1).VAR   = Tnod ;
        MESH(imesh).GAUSS_VECTOR(1).COMP = {'t-x','t-y','t-z'} ;
        MESH(imesh).GAUSS_VECTOR(1).ELEMENTS = ElemLatSurf' + nelem_TOT  ;
        
        
    end
    
end

ndim = size(COOR,2) ;
DATALOC = [] ; 
DATAIN = DefaultField(DATAIN,'TIMELOC',1) ; 
DATALOC.TIMELOC = DATAIN.TIMELOC; 


GidResults2DFE_multi(NameFile_res,ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATALOC);

% cddd = cd ;
% NAMEFILEOPEN =  [cddd,filesep,NameFile_res] ;
% NAMEMESH = [cddd,filesep,NameFile_msh] ;

%  if DATA.COMPRESS_GID == 1
%     disp('Compressing GID FILE ...')
%     TEXT = ['gid -PostResultsToBinary  ',NAMEFILEOPEN, ' ',NAMEFILEOPEN] ;
%     unix(TEXT) ;
%     delete(NameFile_msh)
%     disp('Done')
% end

DATAIN = DefaultField(DATAIN,'NameFile_res',[]) ; 

if isempty(DATAIN.NameFile_res)
    disp('open GID FILE:')
    clipboard('copy',NameFile_res)
    disp(NameFile_res)
    
    DATARUN = DefaultField(DATARUN,'OPEN_GID_ROMRES',0) ;
    if DATARUN.OPEN_GID_ROMRES ==1
        TTTT = ['gidpost ',NameFile_res] ;
        unix(TTTT);
    end
    
end