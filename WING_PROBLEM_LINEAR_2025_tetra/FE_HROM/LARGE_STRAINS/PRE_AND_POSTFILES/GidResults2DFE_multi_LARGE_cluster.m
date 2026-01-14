function GidResults2DFE_multi_LARGE_cluster(NameFile_res,ndime,NODAL_VECTOR,NODAL_SCALAR,MESH,DATA,DATAinpGID,CLUSTER_SEQUENCE)
% INPUT DATA
%dbstop('5')
if nargin == 0
    load('tmp.mat')
elseif nargin == 5
    DATA = [] ; DATAinpGID =[] ;
elseif nargin  == 6
    DATAinpGID =[] ;
end
%%%%%%


fid_res = fopen(NameFile_res,'wt');
fprintf(fid_res,'GiD Post Results File 1.0 \n');
nmesh = length(MESH) ;

for imesh = 1:length(MESH) ;
    MESHLOC = MESH(imesh) ;
    MESHLOC = DefaultField(MESHLOC,'posgp',[]) ;
    MESHLOC = DefaultField(MESHLOC,'TypeElement','Linear') ;
    
    if isempty(MESHLOC.posgp)
        npg = 1;
    else
        npg = size(MESHLOC.posgp,2);
    end
    fprintf(fid_res,['GaussPoints "',MESHLOC.NAMEMESH,'" Elemtype ',MESHLOC.TypeElement,'\n']);
    fprintf(fid_res,['Number of Gauss Points: ',num2str(npg),'\n']);
    fprintf(fid_res,'Nodes not included\n');
    fprintf(fid_res,'Natural Coordinates: Internal\n');
    fprintf(fid_res,'End GaussPoints\n');
    npG(imesh) = npg ;
end


%if ndime == 2
%else
%   nstrain =6 ;
%end
% try
%     DATA = DefaultField(DATA,'NODES',1:nnod) ;
%     nnod = length(DATA.NODES) ;
% catch
%     DATA.NODES = 1:nnod  ;
% end
% stepsSH  =1; % Linear analysis, just one step
% TIMEVECTOR = 1;
DATA = DefaultField(DATA,'TIME_DISCRETIZATION',1) ;
if isempty(DATA.TIME_DISCRETIZATION)
    DATA.TIME_DISCRETIZATION = 1;
end
%DATA = DefaultField(DATA,'STEPS_TO_STORE',1:length(DATA.TIME_DISCRETIZATION)) ;

%DATA = DefaultField(DATA,'STEPS_TO_PRINT', 1:length(DATA.STEPS_TO_STORE)) ;

%DATA.STEPS_TO_PRINT = intersect(DATA.STEPS_TO_PRINT,DATA.STEPS_TO_STORE) ;

DATAinpGID = DefaultField(DATAinpGID,'OPERreconstr',[]) ;
DATAinpGID = DefaultField(DATAinpGID,'ADDITION_VARIABLE',[]) ;

for ivar = 1:length(NODAL_VECTOR)
    DATAinpGID.OPERreconstr = DefaultField(DATAinpGID.OPERreconstr,NODAL_VECTOR(ivar).NAME_SNAP,[]) ;
    DATAinpGID.ADDITION_VARIABLE = DefaultField(DATAinpGID.ADDITION_VARIABLE,NODAL_VECTOR(ivar).NAME_SNAP,[]) ;
end


for istep = 1:length(DATA.TIME_DISCRETIZATION)
    
    time_step = DATA.TIME_DISCRETIZATION(istep) ;
    
    istepGLO = find(time_step == DATA.STEPS) ;
    
    % NODAL VARIABLES
    %%%%%%%%%%%%%%%%%
    % VECTOR
    NAMEFIELDS = [] ;
    for ivar = 1:length(NODAL_VECTOR)
        NAMEFIELDS.VARIABLE = NODAL_VECTOR(ivar).NAME ;
        NAMEFIELDS.COMP = NODAL_VECTOR(ivar).COMP ;
        
        VAR = NODAL_VECTOR(ivar).VAR(:,istep) ;
        if ~isempty( DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP));
            iCL = CLUSTER_SEQUENCE(istepGLO)  ; % Current cluster
            coeff = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).coeff ;
            
            if ~iscell(coeff) &&  coeff == 1
                if  isstruct(DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS)
                    nrowsLOC = size(DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS.coeff{iCL},2) ;
                    VAR = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS.coeff{iCL}*VAR(1:nrowsLOC,:);
                    VAR = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS.BASIS*VAR; 
                else
                    nrowsLOC = size(DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS{iCL},2) ;
                    VAR = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS{iCL}*VAR(1:nrowsLOC,:);
                end
                
            else
                nrowsLOC = size(DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).coeff{iCL},2) ;
                VAR = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).coeff{iCL}*VAR(1:nrowsLOC);
                VAR = DATAinpGID.OPERreconstr.(NODAL_VECTOR(ivar).NAME_SNAP).BASIS*VAR;
            end
            
        end
        
        if ~isempty(DATAinpGID.ADDITION_VARIABLE.(NODAL_VECTOR(ivar).NAME_SNAP))
            VARadd = DATAinpGID.ADDITION_VARIABLE.(NODAL_VECTOR(ivar).NAME_SNAP).a(:,istepGLO) ;
            VARadd =  DATAinpGID.ADDITION_VARIABLE.(NODAL_VECTOR(ivar).NAME_SNAP).U*VARadd ;
            VAR = VAR + VARadd ;
        end
        
        NODES = NODAL_VECTOR(ivar).NODES ;
        GidPostPlotVectNodesNEW(VAR,NODES,time_step,NAMEFIELDS,fid_res,ndime) ;
    end
    % SCALAR
    NAMEFIELDS =  [] ;
    for ivar = 1:length(NODAL_SCALAR)
        NAMEFIELDS.VARIABLE = NODAL_SCALAR(ivar).NAME ;
        VAR = NODAL_SCALAR(ivar).VAR(:,istep) ;
        NODES = NODAL_SCALAR(ivar).NODES' ;
        GidPostPlotScalarNodesNEW(VAR,ndime,NODES,time_step,NAMEFIELDS,fid_res) ;
    end
    
    % GAUSS RESULTS
    % --------------
    for imesh = 1:length(MESH)
        MESHloc =  MESH(imesh) ;  ;
        NAMEMESH_loc = MESHloc.NAMEMESH;
        npg = npG(imesh);
        MESHloc= DefaultField(MESHloc,'GAUSS_SCALAR',[]) ;
        NAMEFIELDS = [] ;
        for ivar = 1:length(MESHloc.GAUSS_SCALAR)
            
            NAMEFIELDS.VARIABLE = MESHloc.GAUSS_SCALAR(ivar).NAME ;
            VAR = MESHloc.GAUSS_SCALAR(ivar).VAR(:,istep) ;
            ELEMENTS = MESHloc.GAUSS_SCALAR(ivar).ELEMENTS ;
            GidPostPlotScalarGaussNEW(VAR,npg,ELEMENTS,time_step,NAMEFIELDS,fid_res,NAMEMESH_loc) ;
        end
        
        MESHloc = DefaultField(MESHloc,'GAUSS_VECTOR',[]) ;
        
        for ivar = 1:length(MESHloc.GAUSS_VECTOR)
            NAMEFIELDS.VARIABLE = MESHloc.GAUSS_VECTOR(ivar).NAME ;
            VAR = MESHloc.GAUSS_VECTOR(ivar).VAR(:,istep) ;
            NAMEFIELDS.COMP = MESHloc.GAUSS_VECTOR(ivar).COMP ;
            ELEMENTS = MESHloc.GAUSS_VECTOR(ivar).ELEMENTS ;
            ndimeLOC = length(NAMEFIELDS.COMP) ;
            ncomp = DATA.MESH.nstrain*npg ;
            VAR = reshape(VAR,ncomp,[]) ;
            GidPostPlotVectorGaussNEW(VAR,ndime,npg,time_step,NAMEFIELDS,fid_res,ndimeLOC,ELEMENTS,NAMEMESH_loc) ;
            
            
            
        end
        MESHloc = DefaultField(MESHloc,'GAUSS_MATRIX',[]) ;
        %     DATAinpGID.OPERreconstr = DefaultField(DATAinpGID.OPERreconstr,NODAL_VECTOR.NAME_SNAP,[]) ;
        for ivar = 1:length(MESHloc.GAUSS_MATRIX)
            NAMEFIELDS.VARIABLE = MESHloc.GAUSS_MATRIX(ivar).NAME ;
            VAR = MESHloc.GAUSS_MATRIX(ivar).VAR(:,istep) ;
            nstrain =  DATA.MESH.nstrainACTUAL ;
            
            ncomp = nstrain*npg ;
            VAR = reshape(VAR,ncomp,[]) ;
            NAMEFIELDS.COMP = MESHloc.GAUSS_MATRIX(ivar).COMP ;
            ELEMENTS = MESHloc.GAUSS_MATRIX(ivar).ELEMENTS ;
            GidPostPlotMatGaussNEWlarge(VAR,nstrain,npg,time_step,NAMEFIELDS,fid_res,ndime,ELEMENTS,NAMEMESH_loc)
            
        end
        
        
        %
        %
        %          MESHLOC.GAUSS_MATRIX(1).NAME = 'STRESS' ;
        % MESHLOC.GAUSS_MATRIX(1).VAR = stressGLO ;
        % MESHLOC.GAUSS_MATRIX(1).COMP = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'}  ;
        % MESHLOC.GAUSS_MATRIX(1).ELEMENTS = IND_ELEM_MESHES{imesh} ;
        %
        
    end
    
    
    
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% STRESSES ........................
    %     % ...................................
    %     if ~isempty(stressGLO)
    %         NAMEFIELDS.VARIABLE = 'STRESS' ;
    %         NAMEFIELDS.COMP  = {'Stress-xx','Stress-yy','Stress-zz','Stress-xy','Stress-yz','Stress-xz'} ;
    %         setElementsRED = [] ;
    %         GidPostPlotMatGauss(stressGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% REDUCED STRESSES ........................
    %     % ...................................
    %     DATA  =DefaultField(DATA,'setElementsRED',[]) ;
    %     if ~isempty(stressGLO) & ~isempty(DATA.setElementsRED);
    %         NAMEFIELDS.VARIABLE = 'GAPPY STRESS' ;
    %         NAMEFIELDS.COMP  = {'StressG-xx','StressG-yy','StressG-zz','StressG-xy','StressG-yz','StressG-xz'} ;
    %         setElementsRED = DATA.setElementsRED ;
    %         var = stressGLO(:,DATA.setElementsRED) ;
    %         GidPostPlotMatGauss(var,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% STRAINS ........................
    %     % ...................................
    %     if ~isempty(strainGLO)
    %         NAMEFIELDS.VARIABLE = 'STRAIN' ;
    %         NAMEFIELDS.COMP  = {'Strain-xx','Strain-yy','Strain-zz','Strain-xy','Strain-yz','Strain-xz'} ;
    %         setElementsRED = [] ;
    %         GidPostPlotMatGauss(strainGLO,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% VON MISES STRESS ........................
    %     % ...................................
    %     DATA = DefaultField(DATA,'PRINT_VONMISES_STRESS',1) ;
    %     DATA = DefaultField(DATA,'stressVONMISES',[]) ;
    %     if ~isempty(DATA.stressVONMISES) & DATA.PRINT_VONMISES_STRESS == 1
    %         NAMEFIELDS.VARIABLE = 'VON_MISES' ;
    %         setElementsRED = [] ;
    %         GidPostPlotScalarGauss(DATA.stressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %     DATA = DefaultField(DATA,'MAXstressVONMISES',[]) ;
    %     if ~isempty(DATA.MAXstressVONMISES)
    %         NAMEFIELDS.VARIABLE = 'MAX_VON_MISES' ;
    %         setElementsRED = [] ;
    %         GidPostPlotScalarGauss(DATA.MAXstressVONMISES,nstrain,npg,nElem,time_step,NAMEFIELDS,DATA,fid_res,ndime,setElementsRED) ;
    %     end
    %
    %     if ~isempty(DATA.Forces1D)
    %
    %         fff = fieldnames(DATA.Forces1D) ;
    %         for ifield = 1:length(fff)
    %             FORCE = DATA.Forces1D.(fff{ifield}) ;
    %             if ~isempty(FORCE)
    %                 NAMEFIELDS.VARIABLE = fff{ifield} ;
    %                 GidPostPlotScalarNodes(FORCE,ndime,nnod,time_step,NAMEFIELDS,DATA,fid_res) ;
    %             end
    %         end
    %     end
    
end




fclose(fid_res);