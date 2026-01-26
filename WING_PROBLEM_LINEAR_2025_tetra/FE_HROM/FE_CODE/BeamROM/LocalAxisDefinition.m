function [RotMatrixGlo Elements2Print] = LocalAxisDefinition(MESH1D,MESH3D,DATAIN)

% ROTATION MATRICES
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp1.mat')
end

[nnode ndim]= size(MESH1D.COOR) ;
[nelem nnodeE]= size(MESH1D.CN) ;
Elements2Print = [] ;
TypeElements2Print = [] ;
SliceElementToPrint = [] ; 
RotMatrixGlo = zeros(ndim,nelem*ndim) ;
for iENTITY  = 1: length(MESH1D.PROP)
    itypeGLO = MESH1D.PROP(iENTITY).INDEX_ENTITY ;
    switch MESH1D.PROP(iENTITY).TYPE ;
        case {'SLICE','BEAM'}
            for jslice = 1:length(itypeGLO)
                itypeSLICE = itypeGLO(jslice) ;
                ROT_LOC_SLICE = MESH3D.SLICES(itypeSLICE).DATA3D.ROTATIONf1f2; % Local rotation slice (Curved elements)
                [RotMatrixGlo Elements2PrintLoc] =RotationMatrixSlices(MESH1D,iENTITY,RotMatrixGlo,ROT_LOC_SLICE,itypeSLICE) ;
                Elements2Print  = [Elements2Print ; Elements2PrintLoc] ;
                SliceElementToPrint = [SliceElementToPrint ; itypeSLICE*ones(size(Elements2PrintLoc))] ; 
                TypeElements2Print = [TypeElements2Print ; iENTITY*ones(size(Elements2PrintLoc))] ;
            end
        case 'JOINT'
            warning('Option not implemented')
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
        RotMatrixGlo,Elements2Print(:,1)') ;
    
    disp(['Open GID file to check direction of local axes: ']) ;
    disp([ cd,filesep,NameFile_res]) ;
end

Elements2Print = [Elements2Print,TypeElements2Print , SliceElementToPrint] ;