function [RotMatrixGlo,ORDER_CONNECTIVITIES,xSIGN] = LocalAxisDefinitionBeams(MESH1D,MESH3D,DATAIN)

% ROTATION MATRICES
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp4.mat')
end

[nnode ndim]= size(MESH1D.COOR) ;
[nelem nnodeE]= size(MESH1D.CN) ;

RotMatrixGlo = zeros(ndim,nelem*ndim) ; % Global matrix containinng the rotation of each SLICE
ORDER_CONNECTIVITIES = ones(nelem,1) ;

for iENTITY  = 1: length(MESH1D.PROP)
    itypeGLO = MESH1D.PROP(iENTITY).INDEX_ENTITY ;
    switch MESH1D.PROP(iENTITY).TYPE ;
        case {'SLICE','BEAM'}
            for jslice = 1:length(itypeGLO)
                itypeSLICE = itypeGLO(jslice) ;
                MESH3D.SLICES(itypeSLICE).DATA3D = DefaultField(MESH3D.SLICES(itypeSLICE).DATA3D,'ROTATIONf1f2',eye(3));
                
                ROT_LOC_SLICE = MESH3D.SLICES(itypeSLICE).DATA3D.ROTATIONf1f2; % Local rotation slice (Curved elements)
              
                
                
                CENTRf1 = MESH3D.SLICES(itypeSLICE).DATA3D.CENTRf1(:); % Centroid f1, reference configuration
                CENTRf2= MESH3D.SLICES(itypeSLICE).DATA3D.CENTRf2(:); % Centroid f2, reference configuration
               
                % Compute the rotation matrix associated to each slice
                [RotMatrixGlo,ORDER_CONNECTIVITIES,MESH1D.CN,xSIGN] =RotationMatrixSlices(MESH1D,iENTITY,...
                    RotMatrixGlo,CENTRf1,CENTRf2,itypeSLICE,ORDER_CONNECTIVITIES,ROT_LOC_SLICE,DATAIN) ;
                
            end
            %         case 'JOINT'
            %             warning('Option not implemented')
    end
end
DATAIN = DefaultField(DATAIN,'PlotNormalsX',0) ;

if DATAIN.PlotNormalsX == 1
    %GID PRINTING
    [dummy1 NAMEMESH_STRUCTURE    dummy2]= fileparts(MESH1D.NAME) ;
    [dummy1 NAMEMESH_STRUCTURE    dummy2]= fileparts(MESH1D.NAME) ;
    NameFile_msh = ['GIDPOST',filesep,NAMEMESH_STRUCTURE,'_',DATAIN.NAME_project,'local_axes','.msh'] ;
    NameFile_res = ['GIDPOST',filesep,NAMEMESH_STRUCTURE,'_',DATAIN.NAME_project,'local_axes','.res'] ;
    
    GidMesh2DFE_multi(NameFile_msh,MESH1D.COOR,{MESH1D.CN},'',{MESH1D.MaterialType},{MESH1D.TypeElement},{'1D'});
    Normals1DScheleton(DATAIN,{MESH1D.TypeElement},{'1D'},NameFile_res,...
        RotMatrixGlo) ;
    
    disp(['Open GID file to check direction of local axes: ']) ;
    disp([ cd,filesep,NameFile_res]) ;
end

