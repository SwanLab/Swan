function [Basis_tPRESS,S_tPRESS] = tPRESS_fromBasisFromDisplacements(BasisU,CASES,NAMEsnap_base,DOFr,DOFl,DATAparamSTUDY,...
    OPERFE,MATPRO,DATAoffline,DATA)

% Computation of follower loads as a function of projected displacements
% JAHO, 5-July-2021, UPCT, Cargagena 


if nargin == 0
    load('tmp.mat')
end


tPRESS_error = zeros(length(CASES),1) ;

SNAP_tPRESS_proj = cell(1,length(CASES)) ;


DATAparamSTUDY = DefaultField(DATAparamSTUDY,'DisableStressComparison',0) ; % Because of the size of the matrix ---too demanding

if DATAparamSTUDY.DisableStressComparison == 1
    StoreBlockSnapshotsForSVD = 1;
else
    StoreBlockSnapshotsForSVD = 0 ;
end

DATAparamSTUDY = DefaultField(DATAparamSTUDY,'StoreBlockSnapshotsForSVD',StoreBlockSnapshotsForSVD) ; % Because keeping it in RAM is too costly


 


for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    iniCLUS  = 1;
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAP_tPRESS_proj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAP_tPRESS_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    disp(['Computing PK2 stresses'])
    disp('Loop over clusters each project')
    for iloc = 1:length(NAME_SNAP_loc)
        disp(['iclus=',num2str(iloc)])
        Nameloc = NAME_SNAP_loc{iloc} ;
        
        if  DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
            [filepath,name,ext] = fileparts(Nameloc) ;
            NAMELOC_tpress = [filepath,filesep,name,'_tpress',ext] ;
           
        end
        
        
        if exist(Nameloc)
            load(Nameloc,'SNAP_cluster') ;
            % ---------------------------------------------------------------------------------------------
            
            if DATAparamSTUDY.DisableStressComparison == 0
                SNAP_tPRESS_LOC{iloc} = bsxfun(@times,SNAP_cluster.tPRESS.U',SNAP_cluster.tPRESS.S)' ;
                SNAP_tPRESS_LOC{iloc} = SNAP_tPRESS_LOC{iloc}*SNAP_cluster.tPRESS.V' ;
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
         %   VAR.DISP = d;  
            
            
            % Projected velocities 
            % ---------------------
            coeff = BasisU'*SNAP_cluster.VEL.U ;
            coeff =bsxfun(@times,coeff',SNAP_cluster.VEL.S)' ;
            coeff = coeff*SNAP_cluster.VEL.V' ;
            
            vL = BasisU*coeff; % Snapshot displacements (local)
            vR = bsxfun(@times,SNAP_cluster.VEL.U(DOFr,:)',SNAP_cluster.VEL.S)' ;
            vR = vR*SNAP_cluster.VEL.V' ;
            ndof = size(vL,1)+size(vR,1) ;
            v = zeros(ndof,size(dL,2)) ;
            v(DOFl,:)  = vL ;
            v(DOFr,:)  = vR ;
          %  VAR.VEL = v;  
            
            %   end
            disp(['Computing tPRESS from displacements ...'])
            tPRESSst = zeros(DATA.MESH.HYDRO.ngausT*DATA.MESH.ndim,size(d,2)) ; 
            for  itime = 1:size(d,2)
               
                VAR.DISP = d(:,itime) ; 
                VAR.VEL = v(:,itime) ; 
                if itime >2 
                DATA.DeltaT  = DATA.STEPS(itime) - DATA.STEPS(itime-1) ; 
                else
                    DATA.DeltaT = DATA.STEPS(itime+1) - DATA.STEPS(itime);
                end
                DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO = 1 ; 
            [~, tPRESSst(:,itime),~ ]= PressureForcesFromDisplacements(DATA,OPERFE,VAR)   ; 
            end
            
            disp(['...Done'])
            
            if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
                % SNAP_tPRESS_proj_LOC{iloc} will be a path to the matrix
                % rather than the matrix itself
                % ----------------------------------------------------------
                %                 NAMELOCPK2 = [filepath,filesep,name,'_PK2',ext] ;
                %                 NAMELOCPK1= [filepath,filesep,name,'_PK1',ext] ;
                save(NAMELOC_tpress,'tPRESSst')  ;
               
                SNAP_tPRESS_proj_LOC{iloc} =NAMELOC_tpress ;
                
            else
                SNAP_tPRESS_proj_LOC{iloc} = tPRESSst  ;
              
            end
        end
    end
    
    if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 0
        
        %else
        SNAP_tPRESS_proj_LOC = cell2mat(SNAP_tPRESS_proj_LOC) ;   % Approximate
    end
    
    
    
    if DATAparamSTUDY.DisableStressComparison == 0
        SNAP_tPRESS_LOC = cell2mat(SNAP_tPRESS_LOC) ;   % exact
        
        % Check if it meets the error criterion
        TIME_STEPS_compare = 2:size(SNAP_tPRESS_LOC,2) ;  % The initial step is excluded
        % Check when stresses are maximum
        % --------------------------------
%         ncompare = 10 ;  % Number of maximum stresses
%         ncompare = min(ncompare,length(TIME_STEPS_compare)) ;
      %  nSTRESSfe  = sqrt(sum(SNAP_tPRESS_LOC(:,TIME_STEPS_compare).^2,1)) ;
      %  [sortSTRESS,INDsorted ] = sort(nSTRESSfe,'descend')   ; % timeMAX is the time at which differences are maximum
      %  Times_compare_max = TIME_STEPS_compare(INDsorted(1:ncompare)) ; % Time steps to compare
      %  DIFF_stresses_atMAX = SNAP_tPRESS_LOC(:,Times_compare_max)-SNAP_tPRESS_proj_LOC(:,Times_compare_max) ;
      %  nDIFF_stresses_atMAX  =norm(DIFF_stresses_atMAX,'fro') ;
      %  nstressATMAX = norm(SNAP_tPRESS_LOC(:,Times_compare_max),'fro') ;
      %  relERRORmax = nDIFF_stresses_atMAX/nstressATMAX ;
        averagerror = norm(SNAP_tPRESS_LOC(:,TIME_STEPS_compare)-SNAP_tPRESS_proj_LOC(:,TIME_STEPS_compare),'fro')/norm(SNAP_tPRESS_LOC(:,TIME_STEPS_compare),'fro')  ;
        tPRESS_error(iproj) = averagerror ;
        disp(['Project = ',num2str(iproj),'; average ERROR tPRESS= ',num2str(averagerror)]);
        if tPRESS_error(iproj) > DATAoffline.errorSTRESS
            %  dbstop('129')
            error('STRESS ERROR CRITERION NOT MET: TRY WITH A LOWER ERROR TOLERANCE FOR DISPLACEMENTS ')
        end
    end
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    %   [UU,SS,VV] =  RSVDT(SNAP_tPRESS_proj_LOC) ;
    SNAP_tPRESS_proj{iproj} = SNAP_tPRESS_proj_LOC ;
    clearvars SNAP_tPRESS_proj_LOC
    % PK1 stresses
    % ----------------
     
    
end

% BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK2 STRESSES (FOR RECONSTRUCTION PURPOSES)
% -------------------------------------------------------------------------------------------------



if  DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
    SNAP_tPRESS_proj =   horzcat( SNAP_tPRESS_proj{:} );
    TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAP_tPRESS_proj),1)' ;
else
    TOL_BLOCK = DATAoffline.errorPK2stress_basis*ones(length(SNAP_tPRESS_proj),1)' ;
end


aaa = cellfun(@isempty,SNAP_tPRESS_proj)  ;
inonempty = find(aaa==0) ; 
SNAP_tPRESS_proj = SNAP_tPRESS_proj(inonempty) ; 

TOL_BLOCK = TOL_BLOCK(inonempty) ; 

disp(['SVD tPRESS snapshots...'])
DATAsvd=[];
DATAsvd.HIDE_OUTPUT = 0 ;
[Basis_tPRESS,S_tPRESS,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAP_tPRESS_proj,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of tPRESS  modes =',num2str(size(Basis_tPRESS,2))])
disp('***********************************************************')


% % BASIS MATRIX FOR THE COLUMN SPACE OF THE SNAPSHOT OF PK1 STRESSES (FOR applying the ECM )
% % ----------------------------------------------------------------------------------------------
% save(DATAoffline.NAMEOFFLINEstore,'BasisStwo') ;
% nmodesS=  size(BasisStwo,2) ;
% clearvars 'BasisStwo'

if DATAparamSTUDY.StoreBlockSnapshotsForSVD == 1
    for iii = 1:length(SNAP_tPRESS_proj)
        delete(SNAP_tPRESS_proj{iii}) ;
    end
end


 
% SNAPstressPonePROJ = SNAPstressPonePROJ(inonempty) ; 
% 
% DATAsvd=[];
% DATAsvd.HIDE_OUTPUT = 0 ;
% [BasisPone,Spone,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(SNAPstressPonePROJ,TOL_BLOCK,DATAsvd) ;
% % disp('***********************************************************')
% % disp(['Number of PK1 stress modes =',num2str(size(BasisPone,2))])
% % disp('***********************************************************')
% % The basis for the PK1 stresses is weighted with the singular values
% 

if ~isempty(S_tPRESS)
Basis_tPRESS = bsxfun(@times,Basis_tPRESS',S_tPRESS)' ;
else
    Basis_tPRESS = [] ; 
end


 