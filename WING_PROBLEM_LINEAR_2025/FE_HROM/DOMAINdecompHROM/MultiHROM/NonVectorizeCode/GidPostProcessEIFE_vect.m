function [NameFile_res ]= GidPostProcessEIFE_vect(COOR,CN,TypeElement,d, STRESSES_FINE,...
    posgp,NameFile_res,MaterialType,DATAIN,DATA,VON_mises_COARSE,IND_ELEM_MESHES,InternalVarStrain,...
    INTERNAL_VAR_FINE,VON_MISES_FINE);
% Post-processing of results using GID (EIFE method)
% See
%dbstop('5')
if nargin==0
    load('tmp.mat')
end


% Writing results file
% -----------------------------------------------------
% Nodal vectors
% -----------------------------------------------------
NODAL_VECTOR(1).NAME = 'DISPLACEMENTS' ;
if  size(COOR,2) == 2
    NODAL_VECTOR(1).COMP = {'X-DISP','Y-DISP'};
    
else
    NODAL_VECTOR(1).COMP = {'X-DISP','Y-DISP','Z-DISP'};
end
NODAL_VECTOR(1).VAR = d ;
NODAL_VECTOR(1).NODES = 1:size(COOR,1) ;

DATA = DefaultField(DATA,'REACTIONS3D',[]) ;

if ~isempty(DATA.REACTIONS3D)
    NODAL_VECTOR(2).NAME = 'REACTIONS-SE' ;
    NODAL_VECTOR(2).COMP = {'X-REAC','Y-REAC','Z-REAC'};
    NODAL_VECTOR(2).VAR = DATA.REACTIONS3D ;
    NODAL_VECTOR(2).NODES = 1:size(COOR,1) ;
end



% VON-MISES STRESSES, coarse
% -------------------------------------
MESH = [] ;
imesh = 1;
ires = 1;
MESH.TypeElement = TypeElement{imesh} ;
MESH.posgp = [] ;
MESH.NAMEMESH = DATA.NAMEMESH{imesh};
MESH.GAUSS_SCALAR(ires).NAME = 'Max.VonMises_coar' ;
MESH.GAUSS_SCALAR(ires).VAR = VON_mises_COARSE ;
MESH.GAUSS_SCALAR(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;

% Internal variable
% -------------------------------------
if ~isempty(InternalVarStrain)
ires = 2;
MESH.GAUSS_SCALAR(ires).NAME = 'Max.Int.Var.Strain' ;
MESH.GAUSS_SCALAR(ires).VAR = InternalVarStrain  ; 
MESH.GAUSS_SCALAR(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;
end



% if size(stressGLO_coarse.AVERAGE,1) == 6
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
% elseif   size(stressGLO_coarse.AVERAGE,1) == 4
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy'}  ;
% else
%     error('Option not implemented')
% end


% 
% % COARSE-SCALE STRESSES
% % -------------------------------------
% MESH = [] ;
% imesh = 1;
% ires = 1;
% MESH.TypeElement = TypeElement{imesh} ;
% MESH.posgp = posgp{imesh} ;
% MESH.NAMEMESH = DATA.NAMEMESH{imesh};
% MESH.GAUSS_MATRIX(ires).NAME = 'AVG.STRESS' ;
% MESH.GAUSS_MATRIX(ires).VAR = stressGLO_coarse.AVERAGE
% MESH.GAUSS_MATRIX(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;
% if size(stressGLO_coarse.AVERAGE,1) == 6
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
% elseif   size(stressGLO_coarse.AVERAGE,1) == 4
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy'}  ;
% else
%     error('Option not implemented')
% end
% 
% ires = 2;
% MESH.GAUSS_MATRIX(ires).NAME = 'MAX.STRESS' ;
% MESH.GAUSS_MATRIX(ires).VAR = stressGLO_coarse.MAXIM ;
% MESH.GAUSS_MATRIX(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;
% if size(stressGLO_coarse.MAXIM,1) == 6
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
% elseif   size(stressGLO_coarse.MAXIM,1) == 4
%     MESH(imesh).GAUSS_MATRIX(ires).COMP = {'Stress-xx','Stress-yy','Stress-xy','Stress-zz'}  ;
% else
%     error('Option not implemented')
% end









% ------------------------------------------------------
%%%%%%%% VARIABLES DEFINED ON 3D GAUSS POINT MESHES
% -----------------------------------------------------

%MESH = [] ;
NODAL_SCALAR = [] ;
imesh1D_add = 1;
for imeshLOC  = 1:length(STRESSES_FINE)
    imesh = imeshLOC +imesh1D_add;   % The first mesh is 1D
    
    if ~isempty(CN{imesh})
        % ----------
        % 3D mesh
        % ----------
        MESH(imesh).TypeElement = TypeElement{imesh} ;
        MESH(imesh).posgp = posgp{imesh} ;
        MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
        %  nelemE(imesh) = size(CN{imesh},2) ;
        
        %   iRESULTS = imesh -imesh1D_add;
        
        %         if ~isempty(STRESSDATA.VONMISES)
        %             MESH(imesh).GAUSS_SCALAR(1).NAME = 'Von_Mises_3D' ;
        %             MESH(imesh).GAUSS_SCALAR(1).VAR = STRESSDATA.VONMISES{iRESULTS}' ;
        %             MESH(imesh).GAUSS_SCALAR(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
        %         end
        %
        
        PLOT_STRESS_FINE = 0 ; 
        
        if ~isempty(STRESSES_FINE{imeshLOC}) &&  PLOT_STRESS_FINE == 1
            MESH(imesh).GAUSS_MATRIX(1).NAME = 'STRESSfine' ;
            MESH(imesh).GAUSS_MATRIX(1).VAR = STRESSES_FINE{imeshLOC} ;
            MESH(imesh).GAUSS_MATRIX(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
            if DATA.MESH.nstrain == 6
                MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
                
            elseif  DATA.MESH.nstrain == 4
                MESH(imesh).GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-xy','Stress-zz'}  ;
            else
                error('Option not implemented')
            end
            
        end
        
      ires = 0 ; 
        if ~isempty(VON_MISES_FINE{imeshLOC})
            ires = ires+1 ; 
        MESH(imesh).GAUSS_SCALAR(ires).NAME = 'Von Mises FINE' ;
        MESH(imesh).GAUSS_SCALAR(ires).VAR = VON_MISES_FINE{imeshLOC}  ;
        MESH(imesh).GAUSS_SCALAR(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;
        end
        
       
        if ~isempty(INTERNAL_VAR_FINE{imeshLOC})
               ires  = ires+1; 
        MESH(imesh).GAUSS_SCALAR(ires).NAME = 'Int.Var.Str. FINE' ;
        MESH(imesh).GAUSS_SCALAR(ires).VAR = INTERNAL_VAR_FINE{imeshLOC}  ;
        MESH(imesh).GAUSS_SCALAR(ires).ELEMENTS = IND_ELEM_MESHES{imesh} ;
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
%         [Tnod, ElemLatSurf]= TractionForcesPrint(DATA_REFMESH,FORCES_2_PRINT,DATAIN) ;
%         
%         % External traction forces
%         % ----------------- %
%         MESH(imesh).GAUSS_VECTOR(1).NAME = 'TRACTION' ;
%         MESH(imesh).GAUSS_VECTOR(1).VAR   = Tnod ;
%         MESH(imesh).GAUSS_VECTOR(1).COMP = {'t-x','t-y','t-z'} ;
%         MESH(imesh).GAUSS_VECTOR(1).ELEMENTS = ElemLatSurf' + nelem_TOT  ;
%         
        
    end
    
end

ndim = size(COOR,2) ;
% DATALOC = [] ;
% DATAIN = DefaultField(DATAIN,'TIMELOC',1) ;
DATALOC.TIMELOC = DATA.TIME_PRINT;


GidResults2DFE_EIFE_vect(NameFile_res,ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATALOC);

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
    %  clipboard('copy',NameFile_res)
    disp(NameFile_res)
    
      DATAIN = DefaultField(DATAIN,'OPEN_GID_POST',0) ;

         if DATAIN.OPEN_GID_POST ==1
             TTTT = ['gidpost ',NameFile_res] ;
             unix(TTTT);
         end
    
end