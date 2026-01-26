function [BasisPone,Spone,nmodesS] = StressesBasisFromDisplacements(BasisU,CASES,NAMEsnap_base,DOFr,DOFl,DATAparamSTUDY,...
    OPERFE,MATPRO,DATAoffline,DATA)

% Computation of BasisPone (PK1 stress basis matrix) and BasisStwo (PK2
% stress basis matrix) corresponding to the displacements projected onto
% the reduced basis matrix
% JAHO, 9-June-2021


if nargin == 0
    load('tmp1.mat')
end


STRESS_PK2_error = zeros(length(CASES),1) ;

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressPonePROJ = cell(1,length(CASES)) ;


DATAparamSTUDY = DefaultField(DATAparamSTUDY,'DisableStressComparison',0) ; % Because of the size of the matrix ---too demanding

if DATAparamSTUDY.DisableStressComparison == 1
    StoreBlockSnapshotsForSVD = 1;
else
    StoreBlockSnapshotsForSVD = 0 ;
end

DATAparamSTUDY = DefaultField(DATAparamSTUDY,'StoreBlockSnapshotsForSVD',StoreBlockSnapshotsForSVD) ; % Because keeping it in RAM is too costly


%  if  DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
%
%
%  end
%


for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    iniCLUS  = 1;
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWOproj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    SNAPstressPonePROJ_LOC = cell(1,length(NAME_SNAP_loc)) ;
    disp(['Computing PK2 stresses'])
    disp('Loop over clusters each project')
    for iloc = 1:length(NAME_SNAP_loc)
        disp(['iclus=',num2str(iloc)])
        Nameloc = NAME_SNAP_loc{iloc} ;
        
        if  DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
            [filepath,name,ext] = fileparts(Nameloc) ;
            NAMELOCPK2 = [filepath,filesep,name,'_PK2',ext] ;
            NAMELOCPK1= [filepath,filesep,name,'_PK1',ext] ;
        end
        
        
        if exist(Nameloc)
            load(Nameloc,'SNAP_cluster') ;
            % ---------------------------------------------------------------------------------------------
            
            if DATAparamSTUDY.DisableStressComparison == 0
                SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U',SNAP_cluster.PK2STRESS.S)' ;
                SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
            end
            
            
            % Projected displacements
            % -----------------------
            coeff = BasisU'*SNAP_cluster.DISP.U(DOFl,:) ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
            coeff = coeff*SNAP_cluster.DISP.V' ;
            % ------------------------------------------------
            dL = BasisU*coeff; % Snapshot displacements (local)
            dR = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
            dR = dR*SNAP_cluster.DISP.V' ;
            ndof = size(dL,1)+size(dR,1) ;
            d = zeros(ndof,size(dL,2)) ;
            d(DOFl,:)  = dL ;
            d(DOFr,:)  = dR ;
            %   end
            disp(['Computing stresses from displacements ...'])
            [~,STRESSPK2,~,~,STRESSPK1,~] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;
            disp(['...Done'])
            
            if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
                % SNAPstressSTWOproj_LOC{iloc} will be a path to the matrix
                % rather than the matrix itself
                % ----------------------------------------------------------
                %                 NAMELOCPK2 = [filepath,filesep,name,'_PK2',ext] ;
                %                 NAMELOCPK1= [filepath,filesep,name,'_PK1',ext] ;
                save(NAMELOCPK2,'STRESSPK2')  ;
                save(NAMELOCPK1,'STRESSPK1')  ;
                SNAPstressSTWOproj_LOC{iloc} =NAMELOCPK2 ;
                SNAPstressPonePROJ_LOC{iloc} = NAMELOCPK1 ;
            else
                SNAPstressSTWOproj_LOC{iloc} = STRESSPK2  ;
                SNAPstressPonePROJ_LOC{iloc} = STRESSPK1  ;
            end
        end
    end
    
    if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 0
        
        %else
        SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    end
    
    
    
    if DATAparamSTUDY.DisableStressComparison == 0
        SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact
        
        % Check if it meets the error criterion
        TIME_STEPS_compare = 2:size(SNAPstressSTWO_LOC,2) ;  % The initial step is excluded
        % Check when stresses are maximum
        % --------------------------------
        ncompare = 10 ;  % Number of maximum stresses
        ncompare = min(ncompare,length(TIME_STEPS_compare)) ;
        nSTRESSfe  = sqrt(sum(SNAPstressSTWO_LOC(:,TIME_STEPS_compare).^2,1)) ;
        [sortSTRESS,INDsorted ] = sort(nSTRESSfe,'descend')   ; % timeMAX is the time at which differences are maximum
        Times_compare_max = TIME_STEPS_compare(INDsorted(1:ncompare)) ; % Time steps to compare
        DIFF_stresses_atMAX = SNAPstressSTWO_LOC(:,Times_compare_max)-SNAPstressSTWOproj_LOC(:,Times_compare_max) ;
        nDIFF_stresses_atMAX  =norm(DIFF_stresses_atMAX,'fro') ;
        nstressATMAX = norm(SNAPstressSTWO_LOC(:,Times_compare_max),'fro') ;
        relERRORmax = nDIFF_stresses_atMAX/nstressATMAX ;
        STRESS_PK2_error_average = norm(SNAPstressSTWO_LOC(:,TIME_STEPS_compare)-SNAPstressSTWOproj_LOC(:,TIME_STEPS_compare),'fro')/norm(SNAPstressSTWO_LOC(:,TIME_STEPS_compare),'fro')  ;
        STRESS_PK2_error(iproj) = relERRORmax ;
        disp(['Project = ',num2str(iproj),'; average ERROR stress PK2= ',num2str(STRESS_PK2_error_average)]);
        disp(['Project = ',num2str(iproj),';  ERROR stress maximum PK2= ',num2str(relERRORmax),' (max at m=',num2str(ncompare),' time steps)']);
        if STRESS_PK2_error(iproj) > DATAoffline.errorSTRESS
            %  dbstop('129')
            error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
        end
    end
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    %   [UU,SS,VV] =  RSVDT(SNAPstressSTWOproj_LOC) ;
    SNAPstressSTWOproj{iproj} = SNAPstressSTWOproj_LOC ;
    clearvars SNAPstressSTWOproj_LOC
    % PK1 stresses
    % ----------------
    if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
        SNAPstressPonePROJ{iproj} = SNAPstressPonePROJ_LOC ;
    else
        SNAPstressPonePROJ{iproj} = cell2mat(SNAPstressPonePROJ_LOC) ;   % exact
    end
    
end

% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% -------------------------------------------------------------------------------------------------



if  DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
    SNAPstressPonePROJ =   horzcat( SNAPstressPonePROJ{:} );
    SNAPstressSTWOproj =   horzcat( SNAPstressSTWOproj{:} );
    TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;
else
    TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAPstressSTWOproj),1)' ;
end


aaa = cellfun(@isempty,SNAPstressSTWOproj)  ;
inonempty = find(aaa==0) ; 
SNAPstressSTWOproj = SNAPstressSTWOproj(inonempty) ; 

TOL_BLOCK = TOL_BLOCK(inonempty) ; 

disp(['SVD stress snapshots...'])
DATAsvd=[];
DATAsvd.HIDE_OUTPUT = 0 ;
[BasisStwo,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressSTWOproj,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of PK2 stress modes =',num2str(size(BasisStwo,2))])
disp('***********************************************************')
% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
% ----------------------------------------------------------------------------------------------
save(DATAoffline.NAMEOFFLINEstore,'BasisStwo') ;
nmodesS=  size(BasisStwo,2) ;
clearvars 'BasisStwo'

if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
    for iii = 1:length(SNAPstressSTWOproj)
        delete(SNAPstressSTWOproj{iii}) ;
    end
end


 
SNAPstressPonePROJ = SNAPstressPonePROJ(inonempty) ; 

DATAsvd=[];
DATAsvd.HIDE_OUTPUT = 0 ;
[BasisPone,Spone,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
% disp('***********************************************************')
% disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
% disp('***********************************************************')
% The basis for the PK1 stresses is weighted with the singular values



BasisPone = bsxfun(@times,BasisPone',Spone)' ;


if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
    for iii = 1:length(SNAPstressPonePROJ)
        delete(SNAPstressPonePROJ{iii}) ;
    end
end