function  UfluctFACE =  PrepareFluctuationsEIFEM(INPUT_PROBLEMS,DATAoffline,dBoundary_FLUCT,...
    MESH,MATPRO,OPERFE,DISP_CONDITIONS,OTHER_output) ;
% JAHO, 4-Dec-2023
% Barcelona, Balmes 185.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/14_GIVEN_FLUCT_Q4.mlx
if nargin == 0
    load('tmp.mat')
end

% geometric MASS MATRIX INTERFACE boundary
% Patterned after  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/FaceNOCornersDOFS_Q4.m
DATAcommon = feval(INPUT_PROBLEMS) ;

% Boundary DOFS
ndim  = size(MESH.COOR,2) ;
bDOFS = MESH.INFO_PERIODIC_CONDITIONS.boundaryDOFS ;  % Boundary DOFs
rnodLOC = large2small(bDOFS,ndim);

ndimLOC = 1;
wDIAGb = CompWeightDiag(OTHER_output.GEOproperties.wSTb,ndimLOC)  ;
Mbound = OTHER_output.GEOproperties.NstB_left'*(wDIAGb*OTHER_output.GEOproperties.NstB_left) ;
Mff = Mbound(rnodLOC,rnodLOC) ;
% This is a matrix relating scalar nodal vectors. Extension to vectors of
% ndim entries per node
nnodeB = size(Mff,1) ;
Mbar = sparse(nnodeB*ndim,nnodeB*ndim) ;
for idim = 1:ndim
    Mbar(idim:ndim:end,idim:ndim:end)  = Mff ;
end

Rb = OTHER_output.GEOproperties.RIGID_BODY_MODES(bDOFS,:) ;

% SVD OF dBoundary_FLUCT
DATALOC1 = [] ;
DATALOC1.TOL = 1e-6 ;
DATALOC1.RELATIVE_SVD = 1; % = DefaultField(DATALOC,'RELATIVE_SVD',1) ;
[Ufluct,S,V]  = WSVDT(cell2mat(dBoundary_FLUCT),Mbar,DATALOC1) ;


BoundaryNodes =rnodLOC ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;

CNb = cell(length(MESH.Indexes_faces_bnd_element),1) ;
for iface  = 1:length(MESH.Indexes_faces_bnd_element)
    CNb{iface} = MESH.CNb(MESH.Indexes_faces_bnd_element{iface},:) ;
end
CNb = cell2mat(CNb(:)) ;
CNbREN  =  RenumberConnectivities( CNb,1:length(BoundaryNodes) );


NameLoc =     ['IMPOSED_DISPT_AND_FLUCT_ALL' ] ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_BASE =[DATAcommon.NameParamStudyLOC ];

NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,Ufluct,...
    [],NameFileMesh,NameFile_res,[],DATALOC);


%%%%%%%% DECOMPOSING FLUCTUATIONS FOR EACH FACE

DATAcommon = DefaultField(DATAcommon,'TypeFunctionDisplacementInterfaces','QUADRILATERAL_LINEAR') ;
switch DATAcommon.TypeFunctionDisplacementInterfaces
    case 'QUADRILATERAL_LINEAR'
        DATAcommon = DefaultField(DATAcommon,'PERIODICITY_FACES',{[1,3],[2,4]}) ;
        % FLUCTUATION FACE 1
        UfluctFACE = cell(1,2) ;
        UfluctALL = cell(1,4) ;
        for iface = 1:length(UfluctFACE)
            indLOC = DATAcommon.PERIODICITY_FACES{iface}(1) ;  % Either 1 or 2
            Floc = MESH.INFO_PERIODIC_CONDITIONS.facesDOFS_loc{indLOC};  % Local DOFs (boundary)
            Uloc = Ufluct(Floc,:);
            TOL = 1e-4 ;
            DATALOC2.RELATIVE_SVD = 1 ;
            [UU,SS,VV] = SVDT(Uloc,TOL,DATALOC2) ;
            
            UfluctFACE{iface}  = UU ;
            
            Floc_plot = MESH.INFO_PERIODIC_CONDITIONS.facesDOFS_loc_withCORNERS{indLOC};  % Local DOFs (boundary)
            Uloc_plot = Ufluct(Floc_plot,:);
            TOL = 1e-4 ;
            DATALOC2.RELATIVE_SVD = 1 ;
            [UU_plot,SS,VV] = SVDT(Uloc_plot,TOL,DATALOC2) ;
            % For ploting purposes
            nrows  = size(Ufluct,1) ;
            for jfacePER = 1:length(DATAcommon.PERIODICITY_FACES{iface})
                jfacePERglo  = DATAcommon.PERIODICITY_FACES{iface}(jfacePER) ;
                UfluctALL{jfacePERglo} = zeros(nrows,size(UU,2)) ;
                FLOCj = MESH.INFO_PERIODIC_CONDITIONS.facesDOFS_loc_withCORNERS{jfacePERglo};
                UfluctALL{jfacePERglo}(FLOCj,:) = UU_plot ;
            end
        end
        
        
        
        NameLoc =     ['FluctPerFace' ] ;
        NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
        NAME_BASE =[DATAcommon.NameParamStudyLOC ];
        
        NAME_MODES_DISP = [NAME_MODES_FOLDER,NAME_BASE,NameLoc] ;
        NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
        NameFile_res = [NAME_MODES_DISP,'.res'] ;
        DATALOC = [];
        
        GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,cell2mat(UfluctALL),...
            [],NameFileMesh,NameFile_res,[],DATALOC);
        
        
        
    otherwise
        error('Option not implemented yet')
        
        
end
