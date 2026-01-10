function [SNAPstressSTWO,BasisS_allclusters ]= GetSnapshotsStresses(CASES,NAMEsnap_base,DATAoffline)

if nargin == 0
    load('tmp1.mat')
end


SNAPstressSTWO  = cell(1,length(CASES)) ;

DATAoffline = DefaultField(DATAoffline,'CompressPK2stresses_originalFE',0);

if  DATAoffline.CompressPK2stresses_originalFE == 1
    SNAP_LEFT_SING_VAL = cell(size(CASES)) ;
end


for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    if  DATAoffline.CompressPK2stresses_originalFE == 1
        SNAP_LEFT_SING_VAL_LOC = cell(size(SNAPstressSTWO_LOC)) ;
    end
    
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        
        
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
        if  DATAoffline.CompressPK2stresses_originalFE == 1
        SNAP_LEFT_SING_VAL_LOC{iloc} =  SNAPstressSTWO_LOC{iloc} ;
        end
        
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        
        
        
    end
    SNAPstressSTWO{iproj} = cell2mat(SNAPstressSTWO_LOC) ;
    
    if DATAoffline.CompressPK2stresses_originalFE == 1
        SNAP_LEFT_SING_VAL{iproj} = cell2mat(SNAP_LEFT_SING_VAL_LOC) ;
    end
    
end

SNAPstressSTWO = cell2mat(SNAPstressSTWO) ;


if  DATAoffline.CompressPK2stresses_originalFE == 1
    
    TOL_BLOCK = DATAoffline.TOL_SVD_MATRIX_ALL*ones(length(SNAP_LEFT_SING_VAL),1)' ;
    
    
    DATAsvd=[];
    DATAsvd.HIDE_OUTPUT =1 ;
    
    disp('Computing basis matrix PK2 stresses (by joining all clusters) ')
    
    [BasisS_allclusters,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAP_LEFT_SING_VAL,TOL_BLOCK,DATAsvd) ;
    disp('***********************************************************')
    disp(['Number of PK2 stress modes (global) =',num2str(size(BasisS_allclusters,2))])
    disp('***********************************************************')
    disp(['Singular Values =',num2str(S')]) ; 
    
else
    BasisS_allclusters = [] ; 
    
    
end
