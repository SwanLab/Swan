function [max_disp,STRAINenergy,max_disp_hist,STRAINenergy_hist,DATA_STEPS] = ExtractSTRAINcoarseBUB(BOUNDARY_CONDIDITONS,DATALOC)
if nargin == 0
    load('tmp1.mat')
end
% PLOT REACTIONS


if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% SNAPSHOTS INFO
% This is the function used for producing the training data

DATAparamSTUDY = [] ;
%DATAparamSTUDY.DOF_TO_STUDY_DISPLACEMENT = 1;


[BoundaryConditions ]= feval(BOUNDARY_CONDIDITONS) ; % Trajectories
[CASES ]= 1:length(BoundaryConditions) ;   % ...Using all training trajectories

% %
% DATAoffline.USE_ELEMENT_BASED_APPROACH = 0;
% DATAoffline.errorDISP = 1e-5 ;% 0.5e-6;   % For each block
% DATAoffline.errorSTRESS = 1e-3;   % For each block. Just for check that the basis matrix for displacements is representative
% DATAoffline.errorFINT = 1e-5; %1e-5;% 1e-5; % 1e-5; %1e-4;
% DATAoffline.errorECM = 0 ;
% DATAoffline.errorPK2stress_basis = 0;
%DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints = 1;
%
%
% PLOT_MODES_WITH_DOFR = 0;  % Set it to zero to generate the modes for the  DOFl

NAME_BASE = [BoundaryConditions.NameParamStudy,'_param_'] ;

NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
FE_VARIABLES_NAMEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'1_FEoper.mat'];

% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;

STRAINenergy =cell(1,length(CASES)) ;
KINETICenergy =cell(1,length(CASES)) ;



load(FE_VARIABLES_NAMEstore,'MATPRO','MESH','DISP_CONDITIONS',...
    'OTHER_output') ;

for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if iproj == 1
        
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        DOFl = DISP_CONDITIONS.DOFl ;
        DOFr = DISP_CONDITIONS.DOFr ;
        DATA.FE_VARIABLES_NAMEstore = FE_VARIABLES_NAMEstore ;
    end
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    STRAINenergyLOC = cell(1,length(NAME_SNAP_loc)) ;
    KINETICenergyLOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPreactions =cell(1,length(NAME_SNAP_loc)) ;
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        if exist(Nameloc)
            load(Nameloc,'SNAP_cluster') ;
            % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
            
            %DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
            
            % Or the whole matrix ....
            DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
            DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
            %    REACTIONS_LOC = cell(1,length(NAME_SNAP_loc)) ;
            
            
            if isfield(SNAP_cluster,'STRAIN_ENERGY')
                
                STRAINenergyLOC{iloc} = bsxfun(@times,SNAP_cluster.STRAIN_ENERGY.U',SNAP_cluster.STRAIN_ENERGY.S)' ;
                STRAINenergyLOC{iloc} =  STRAINenergyLOC{iloc}*SNAP_cluster.STRAIN_ENERGY.V' ;
                
                
            end
            
            if isfield(SNAP_cluster,'KINETICenergyLOC')
                
                %                 STRAINenergyLOC{iloc} = bsxfun(@times,SNAP_cluster.STRAIN_ENERGY.U',SNAP_cluster.STRAIN_ENERGY.S)' ;
                %                 STRAINenergyLOC{iloc} =  STRAINenergyLOC{iloc}*SNAP_cluster.STRAIN_ENERGY.V' ;
                
                KINETICenergyLOC{iloc} = bsxfun(@times,SNAP_cluster.KINETIC_ENERGY.U',SNAP_cluster.KINETIC_ENERGY.S)' ;
                KINETICenergyLOC{iloc} =  KINETICenergyLOC{iloc}*SNAP_cluster.KINETIC_ENERGY.V' ;
                ISDYNAMIC = 1 ;
            else
                ISDYNAMIC = 0 ;
            end
            
            
            
            REACTIONS_LOC{iloc} = bsxfun(@times,SNAP_cluster.RESID.U(:,:)',SNAP_cluster.RESID.S)' ;
            REACTIONS_LOC{iloc} =  REACTIONS_LOC{iloc}*SNAP_cluster.RESID.V' ;
            
            
        end
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    KINETICenergy{iproj} = cell2mat(KINETICenergyLOC);
    SNAPreactions{iproj} =  cell2mat(REACTIONS_LOC);
    
    STRAINenergy{iproj} = cell2mat(STRAINenergyLOC);
    
    
    %  SNAPdisp{iproj} =  SNAPdisp{iproj}(DATAparamSTUDY.DOF_TO_STUDY_DISPLACEMENT,:)   ;
    
    
end

SNAPdisp = cell2mat(SNAPdisp) ;


disp(['--------------------------------------------------------'])
disp(['Extracting information coarse-scale strain'])
disp(['--------------------------------------------------------'])
disp(['Element = ',num2str(DATALOC.ElementToStudy)])
disp(['--------------------------------------------------------'])
strainHISTORY = OTHER_output.strainCOARSE_history{DATALOC.ElementToStudy} ; 

figure(39)
hold on

maxSTRAIN = max(max(abs(strainHISTORY))) ; 

title(['Evolution coarse-scale strains, element = ',num2str(DATALOC.ElementToStudy),' (divided by,',num2str(maxSTRAIN),'   )' ])
xlabel(['Step '])
ylabel(['AMplitude (over maximum)'])
h = [] ; 
LL = {} ; 
for imode = 1:size(strainHISTORY,1)
    h(imode) = plot(strainHISTORY(imode,:)/maxSTRAIN) ; 
    LL{imode} = ['Mode =',num2str(imode)] ; 
end

legend(h,LL)
 
 
 save(DATALOC.NAMEWS_save_strainHISTORY,'strainHISTORY') ; 
 


% BUBBLE APPROACH
% WE ARE ONLY INTERESTED IN THE INTERFACE DOFS
% --------------------------------------------
DOFS_bLOC = DATA.MESHextended.DOFS_bLOC;
SNAPdisp = SNAPdisp(DOFS_bLOC,:) ; 

[max_disp,DOF ]= max(abs(SNAPdisp(:,end))) ;

max_disp_hist = SNAPdisp(DOF,:) ;

SNAPreactions = cell2mat(SNAPreactions) ;
SNAPreactions = SNAPreactions(DOFS_bLOC,:) ; 

% WHAT ARE THE NODES OF FACE
%DATALOC.FACE_TO_ANALIZE = 3;
FACE_TO_ANALIZE = DATALOC.FACE_TO_ANALIZE ;
NODES_FACES= MESH.NODES_FACES{FACE_TO_ANALIZE} ;

ndim = size(MESH.COOR,2) ;

DOFs_face = small2large(NODES_FACES,ndim) ;

DOFS_selected = DOFs_face(DATALOC.DOF_TO_ANALIZE:ndim:end) ;

REACTIONS_time = sum(SNAPreactions(DOFS_selected,:),1) ;

[III,JJJ] = find( BoundaryConditions.INPUTS_PARAMETERS ~=0) ;

DISPALL = BoundaryConditions.INPUTS_PARAMETERS;
imposed_disp = linspace(0,DISPALL(III(1)),length(REACTIONS_time)) ;

figure(106)
hold on
xlabel('Imposed displacement (m)')
ylabel('Reaction (MN)')
h = plot(imposed_disp,REACTIONS_time) ; 


DATALOC = DefaultField(DATALOC,'LEGEND_LOC','') ;
legend(h,{DATALOC.LEGEND_LOC})
%
% times  =2 ;
%
% SNAPreactions = SNAPreactions(:,times);
% ndim = size(MESH.COOR,2) ;
% rx = sum(SNAPreactions(1:ndim:end));
% ry = sum(SNAPreactions(2:ndim:end));
%
disp(['RESULTANT_REACT  = ',num2str(REACTIONS_time(end))]) ;
%  disp(['RESULTANT_REACT Y = ',num2str(ry)]) ;
%
%
STRAINenergy   = cell2mat(STRAINenergy) ;
STRAINenergy_hist = STRAINenergy ;
STRAINenergy = STRAINenergy(end) ;
disp(['STRAIN ENERGY  = ',num2str(STRAINenergy)]) ;



figure(206)
hold on
xlabel('Time (s)')
ylabel('Strain energy (MJ)')
h = plot(DATA.STEPS,STRAINenergy_hist) ; 


DATALOC = DefaultField(DATALOC,'LEGEND_LOC','') ;
legend(h,{DATALOC.LEGEND_LOC})


%
%   disp(['MAX DISP  = ',num2str(max_disp)]) ;
% DATA_STEPS = DATA.STEPS;
%
%
%
%
% % DATAparamSTUDY = DefaultField(DATAparamSTUDY,'NAME_WS_LOC_DATA',['DATAWS/RUN1.mat']) ;
% %
% % [aaaaaa] = fileparts(DATAparamSTUDY.NAME_WS_LOC_DATA) ;
% %
% % if  ~exist(aaaaaa)
% %     mkdir(aaaaaa)
% % end
% % TIME_STEPS= DATA.STEPS ;
% %
% % if ISDYNAMIC == 1
% %
% %     SNAPdisp_DOFS = cell2mat(SNAPdisp) ;
% %     KINETICenergy = cell2mat(KINETICenergy) ;
% %     STRAINenergy = cell2mat(STRAINenergy) ;
% %
% %
% %     DOFS_to_study = DATAparamSTUDY.DOF_TO_STUDY_DISPLACEMENT ;
% %
% %
% %     % REACTIONS ALONG TIME --- AXIAL, SHEAR AND MOMENT
% %     % --------------------------------------------------
% %
% %
% %
% %
% %     save(DATAparamSTUDY.NAME_WS_LOC_DATA,'SNAPdisp_DOFS','DOFS_to_study','TIME_STEPS','KINETICenergy','STRAINenergy') ;
% %
% % else
% %     DOFS_to_study = DATAparamSTUDY.DOF_TO_STUDY_DISPLACEMENT ;
% %
% %     SNAPdisp_DOFS = cell2mat(SNAPdisp) ;
% %     DATAparamSTUDY = DefaultField(DATAparamSTUDY,'NAME_WS_LOC_STATIC','DATAWS/STATIC.mat') ;
% %
% %     save(DATAparamSTUDY.NAME_WS_LOC_DATA,'SNAPdisp_DOFS','DOFS_to_study') ;
% %
% % end
% % DATAINPUTloc = [] ;
% % DATAINPUTloc = DefaultField(DATAINPUTloc,'TIME_STEPS',TIME_STEPS) ;
% % DATAINPUTloc = DefaultField(DATAINPUTloc,'ROTATION_ANGLE_time',@(t) (0*t)) ;
% % DATAINPUTloc = DefaultField(DATAINPUTloc,'DATA_stop',[]) ;
% %
% % if ~isempty(DATAINPUTloc.DATA_stop)
% %     % Define ROTATION_ANGLE_time in a piece-wise fashion
% %     DATAINPUTloc.AngleVersusTime_vector =   AngleVersusTime_picewise(DATAINPUTloc)  ;
% % end
% %
% %
% %
% %
% % %
% % %
% % % [RESULTANTS_FE,ANGLES_FE,WEIGHT,MOMENT_WEIGHT ]= ReactionsResultant_rotating_body(cell2mat(SNAPreactions),DOFr,MESH,DATAINPUTloc,OTHER_output.Fbody_U) ;
% % %
% % %
% % % save(DATAparamSTUDY.NAME_WS_LOC_DATA,'WEIGHT','MOMENT_WEIGHT','RESULTANTS_FE','ANGLES_FE','-append') ;
