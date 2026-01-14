% clc
% clear all
% load('tmp1.mat')
% Prescribing Tractions   Boundary Conditions
%  JAHO.
%20-June-2019  ,
% -----------------------------------------------------------------------------

% FIRST STEP ---> Lines and surfaces created with GID
% ---------------------------------------------------
% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NameFileMesh) ;

% We can have either distributed loads or resultant forces (in turn, distributed over the corresponding line/surface)

%    FORCES.DISTRIBUTED.LINE{ilines} = zeros(ndir,1)  ;  % N/m^2 . Uniform along lines
%     FORCES.RESULTANT.LINE{ilines} = zeros(ndir,1)  ;
%
INPUTS_LOC  = DefaultField(INPUTS_LOC,'FORCES',[]) ;
INPUTS_LOC.FORCES  = DefaultField(INPUTS_LOC.FORCES,'DISTRIBUTED',[]) ;
INPUTS_LOC.RESULTANT  = DefaultField(INPUTS_LOC.FORCES,'RESULTANT',[]) ;
INPUTS_LOC.FORCES.DISTRIBUTED  = DefaultField(INPUTS_LOC.FORCES.DISTRIBUTED,'LINE',[]) ;
INPUTS_LOC.FORCES.DISTRIBUTED  = DefaultField(INPUTS_LOC.FORCES.DISTRIBUTED,'FACE',[]) ;
INPUTS_LOC.FORCES.RESULTANT  = DefaultField(INPUTS_LOC.FORCES.RESULTANT,'LINE',[]) ;
INPUTS_LOC.FORCES.RESULTANT  = DefaultField(INPUTS_LOC.FORCES.RESULTANT,'FACE',[]) ;


if size(COOR,2)==2
    NODES_ENTITIES = NODES_LINES ;
    nnn = min([length(NODES_ENTITIES),length(INPUTS_LOC.FORCES.RESULTANT.LINE)]) ;
    LABEL = 'LINE' ;
else
    NODES_ENTITIES = NODES_FACES ;
    nnn = min([length(NODES_ENTITIES),length(INPUTS_LOC.FORCES.RESULTANT.FACE)]) ;
    LABEL = 'FACE' ;
end

if ~isempty(INPUTS_LOC.FORCES.RESULTANT.(LABEL))
    FORCE_LOC = INPUTS_LOC.FORCES.RESULTANT.(LABEL) ;
    [RRR,CCC] = cellfun(@size,FORCE_LOC);
    Fpnt = zeros(prod(size(COOR)),1) ;
    for iENT = 1:length(FORCE_LOC)
        Fpnt_loc = FORCE_LOC{iENT} ;
        if  any(Fpnt_loc)
            NODES = NODES_ENTITIES{iENT} ;
            Fpnt_loc =  PointLoadsENTITY(NODES,COOR,CONNECTb,...
                TypeElementB,Fpnt_loc,DATA)  ;
            Fpnt = Fpnt + Fpnt_loc ;
        end
    end
    
    
end

DATA.NODES_ENTITIES = NODES_ENTITIES  ; 
 CNb = cell(1,ndim) ; Tnod= cell(1,ndim) ; 
if ~isempty(INPUTS_LOC.FORCES.DISTRIBUTED.(LABEL))
   
    
    FORCE_LOC = INPUTS_LOC.FORCES.DISTRIBUTED.(LABEL) ;
    [RRR] = cellfun(@any,FORCE_LOC);
    
    if any(RRR)
        
        
        
        
        FORCEtypep = INPUTS_LOC.FORCES.DISTRIBUTED ;
        FORCE_LOC = INPUTS_LOC.FORCES.DISTRIBUTED.(LABEL) ; % LABEL = LINE or FACE
        FORCEtypep = DefaultField(FORCEtypep,'ISLOCAL',zeros(size(FORCE_LOC))) ;
        ISLOCAL = FORCEtypep.ISLOCAL ;
        [RRR,CCC] = cellfun(@size,FORCE_LOC);  % Check that there are non-zero forces
        for iENT = 1:length(FORCE_LOC)   % Loop over vorces
            LOAD = FORCE_LOC{iENT} ;
            if  any(LOAD)
                NODES = NODES_ENTITIES{iENT} ;
                [CNbLOC,~] =  ElemBnd(CONNECTb,NODES) ;
                if   ISLOCAL{iENT}
                    % Compute normals to each boundary element
                    NORMALSv = NormalsBoundary(COOR,CNbLOC)  ;      
                     LOAD = LOAD(1)*ones(1,3) ; % To cheat the program
            % The input should be given in the first entry
            % And what about TLOC ? We have to arrange NORMALSv
            % somehow... But how ?
                end
                for idim = 1:length(CNb)
                    tCOMP = LOAD(idim) ;
                    if tCOMP ~= 0
                        CNb{idim} = [CNb{idim} ; CNbLOC] ;
                        if ISLOCAL{iENT} == 0
                            TLOC = tCOMP*ones(size(CNbLOC)) ;
                        else
                            % Local coordinates  (force in the normal direction)
                            % -----------------
                            TLOC = tCOMP*repmat(NORMALSv(idim,:)',1,size(CNbLOC,2)) ;
                            %    TLOC =reshape(TLOC,size(CNbLOC,1),[]) ;
                        end
                        Tnod{idim} = [  Tnod{idim} ; TLOC]  ;
                    end
                end
                
                %                 Fpnt_loc =  PointLoadsENTITY(NODES,COOR,CONNECTb,...
                %                     TypeElementB,Fpnt_loc,DATA)  ;
                %                 Fpnt = Fpnt + Fpnt_loc ;
            end
        end
        
    end
    
end







%
%
% for ilines = 1:nnn
%
%      a = INPUTS_LOC.DISP.LINE{ilines} ;
%     NODES = NODES_ENTITIES{ilines} ;
%     [Gloc,uBARloc,DOFrLOC,DOFmLOC] = ...
%     PRESCRIBED_DISP_OVER_ENTITY(a,NODES,COOR,CONNECTb,TypeElementB,...
%     DATA)  ;
%     G{end+1} =  Gloc ;
%     uBAR{end+1} = uBARloc ;
%     DOFr{end+1} = DOFrLOC ;
%     DOFm{end+1} = DOFmLOC;
%
% end
%
%
% DOFr = cell2mat(DOFr') ;
% uBAR = cell2mat(uBAR') ;
% [DOFrUNIQUE,IND]  = unique(DOFr) ;
%
% if length(DOFrUNIQUE) < length(DOFr)
%
%     DOFm = cell2mat(DOFm') ;
%
%     if ~isempty(DOFm)
%         error('Non-compatible option. Repeated DOFs when imposing BCs cannot exist with affine BCs')
%     end
%
%     DOFr = DOFrUNIQUE;
%     G = [] ;
%
%     uBAR = uBAR(IND) ;
%
% else
%     DOFm = cell2mat(DOFm') ;
%     G = blkdiag(G(:)) ;
%
%
% end
%
%
%
%
% %
% %  MeshLOC_coarse',INPUTS_LOC.NameFileMesh) ;
% % %
% % %
% % %OKK = strcmp(INPUTS_LOC.NameFileMeshLOC_coarse,INPUTS_LOC.NameFileMesh) ;
% % %if OKK ==1
% %     [Gb,dR,DOFr,DOFm] = FIXED_FACES_RVES_fun(a,DOMAINVAR,COOR,CONNECTb,TypeElementB,DATA) ;
% % %else
% %
% % %    error('Option not implemented')
% %
% % %     [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB,DATAcoarsefine] = ...
% % %         PERIODIC_BEAMS_COARSE_FINE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,...
% % %         TypeElementB,INPUTS_LOC.NameFileMeshLOC_coarse,DATA) ;
% %
% % %end
% %
% %
% % %DOMAINVAR.RigidBodyModes =R ;
% %
% %
