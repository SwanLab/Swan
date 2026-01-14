function [stressDOM,reactDOM,Ndom,wDIAG_RHS,Bdom,Wdom,K,CgloDOM] = ...
    ReactionsDOMAIN_new(Nst,fNOD,ListGauss,ndim,DOFs,Fpnt,CNb,NODESrve,Tnod,...
    TypeElementB,CONNECTb,...
    Bst,dLOC,wSTs_RHS,wSTs,Cglo,COOR,nnodeDOM,faceDOFS,BasisUrb,nstrain,idom,iproject,...
    DATAIN,stressGLO,DATA_INPUT_FE,timestepINCLUDE)
% Computing reactions and stresses at each domain
if nargin == 0
    load('tmp2.mat')
end


% STEP 3
%% Computing body forces at domain idom
% -------------------------------------
fbodyDOM = 0 ; Ndom = [] ; wDIAG_RHS = [] ;
if (~isempty(fNOD) && any(fNOD)) || (iproject ==1 && idom ==1)
    % Set of Gauss points associated to domain idom
    ListGaussDOF =  small2large(ListGauss,ndim) ;
    %   if iproject ==1 && idom ==1
    if ~isempty(Nst)
        Ndom = Nst(ListGaussDOF,DOFs) ;
    else
        load(DATAIN.NAME_WS_MODES,'Ndom') ;
    end
    %  end
    WdomN = wSTs_RHS(ListGauss)  ;
    fbodyNODES = fNOD(DOFs) ;
    % Converging WdomN into a diagonal matrix affecting all DOFs
    wDIAG_RHS = CompWeightDiag(WdomN,ndim)  ;
    fbodyDOM = (wDIAG_RHS*Ndom)'*(Ndom*fbodyNODES) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing traction forces -->  FtracDOM
% ------------------------------------------------------
% WE shall invoke function
%  Ftrac = FtracCOMPvect(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb) ;
% but adapting the inputs for domain idom

FpntDOM = zeros(size(Fpnt))  ; % Point loads
FpntDOM(DOFs) = Fpnt(DOFs)  ; % Point loads

CNbDOM = cell(size(CNb)) ;
TnodDOM = cell(size(CNbDOM)) ;
%dbstop('57')


for idim = 1:ndim
    if ~isempty(CNb{idim})
        if  ~isempty(CNb{idim})
            [CNbDOM{idim}, ListEbndDOM]= ElemBnd(CNb{idim},NODESrve) ;
            TnodDOM{idim} = Tnod{idim}(ListEbndDOM,:) ;
        end
    end
end
% OLD_VERSION = 0 ;
% if OLD_VERSION ==1
%     ftracDOM = FtracCOMPvect(COOR,CNbDOM,TypeElementB,FpntDOM,TnodDOM,CONNECTb) ;
%     ftracDOM = ftracDOM(DOFs) ;
% else
% error('Test it !!! ')
DATAIN = DefaultField(DATAIN,'INCLUDE_Fpnt_in_Computing_Traction_Forces',0) ;

if DATAIN.INCLUDE_Fpnt_in_Computing_Traction_Forces == 0
    FpntDOM = zeros(size(FpntDOM)) ;
end



ftracDOM = FtracCOMPvect(COOR,CNbDOM,TypeElementB,FpntDOM,TnodDOM,CONNECTb(idom,:)) ;
ftracDOM = ftracDOM(DOFs) ;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of internal forces
Wdom =  wSTs(ListGauss);

%dbstop('77')
wDIAG = CompWeightDiag(Wdom,nstrain)  ;
ListGaussDOFS =  small2large(ListGauss,nstrain) ;
%if iproject ==1 && idom ==1
if   ~isempty(Bst)
    Bdom = Bst(ListGaussDOFS,DOFs) ;
else
    load(DATAIN.NAME_WS_MODES,'Bdom') ;
end
%else
%   Bdom = Bst ;
%end
%load('tmp0.mat')
%dbstop('82')

%DATAIN = DefaultField(DATAIN,'ComputeInternalForcesFromStresses',0) ;
nsteps = size(dLOC,2) ;

% Computing stresses and internal forces
% ------------------------
CgloDOM = [] ;
if  ~isempty(Cglo)
    CgloDOM=  Cglo(ListGaussDOFS,ListGaussDOFS) ;
else
    load(DATAIN.NAME_WS_MODES,'CgloDOM') ;
end
% if DATAIN.ComputeInternalForcesFromStresses == 0 && nsteps ==1
%     % Stresses arising from displacements (elastic range)
%     strainDOM = Bdom*dLOC ;
%     %strainDOM = Bdom*d(DOFs) ;
%     
%     stressDOM =CgloDOM*strainDOM ; %
%     for istrain =1:nstrain
%         stressDOM(istrain:nstrain:end) = stressDOM(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
%     end
%     %dbstop('90')
% else
    % STresses obtained from the snapshot matrix
    stressDOM = stressGLO(ListGaussDOFS,:) ;
    
%end
fintDOM = (wDIAG*Bdom)'*stressDOM ;



K = [] ;  % Stiffness matrix
if idom == 1 && iproject==1
    K = (Cglo(ListGaussDOFS,ListGaussDOFS)*Bdom) ;
    K = Bdom'*K ;
end

%%%% Reactions
% Traction and body forces in the nonlinear regime (more than one time step)
% This is defined  ExtForces_Pdisp_TIME.m  --- 28-May-2019

DATA_INPUT_FE = DefaultField(DATAIN,'FACTOR_TIME_TRACTION_FORCES',1) ; 
ftracDOM = bsxfun(@times,ftracDOM',DATA_INPUT_FE.FACTOR_TIME_TRACTION_FORCES')' ;
ftracDOM = ftracDOM(:,timestepINCLUDE) ;

if length(fbodyDOM) >1
    DATA_INPUT_FE = DefaultField(DATAIN,'FACTOR_TIME_BODY_FORCES',1) ; 

    fbodyDOM = bsxfun(@times,fbodyDOM',DATA_INPUT_FE.FACTOR_TIME_BODY_FORCES')' ;
    fbodyDOM = fbodyDOM(:,timestepINCLUDE) ;
end
reactDOMorig = fintDOM - (ftracDOM+ fbodyDOM ) ;




%%% Reaction vectors at faces f1 and f2
% f1 = > NODESfaces{idom}{3}
% f2 = > NODESfaces{idom}{1}
%dbstop('100')
%DATAIN.TYPE_REPETITION = 'MULTIDIRECTION';

%-------------------
% 1-D tiled problems
% % ------------------
% f1 = NODESfaces{1} ;
% f2 = NODESfaces{3} ;
% f1 = small2large(f1,ndim) ;
% f2 = small2large(f2,ndim) ;
% f = [f1; f2] ;



% SELF-EQUILIBRATED MODES
% -----------------------
% %%% We seek self-equilibrated reaction force vectors
% Interior nodes are   not accounted for9
DOFtot = 1:ndim*nnodeDOM ;
DOFi = setdiff(DOFtot,faceDOFS) ;
BasisRrb = BasisUrb ;
BasisRrb(DOFi,:) = 0 ;
% Reaction forces with resultant equal to zero
% ---------------------------------------------
%dbstop('145')
reactDOM = zeros(size(reactDOMorig)) ;
reactDOM(faceDOFS) = reactDOMorig(faceDOFS)  - BasisRrb(faceDOFS,:)*(BasisRrb(faceDOFS,:)\reactDOMorig(faceDOFS)) ;



%% ROTATIONS
% -----------
%%%%%
ROTATION =   DATAIN.ROTATION ;  % Maps local to global
if ~isempty(ROTATION)
%     reactDOM = reshape(reactDOM,ndim,[]) ;
%     reactDOM = (ROTATION'*reactDOM);
%     reactDOM  = reactDOM(:) ;
%     
    for istep = 1:size(reactDOM,2)
        reactDOM_istep = reshape(reactDOM(:,istep),ndim,[]) ;
        reactDOM_istep = (ROTATION'*reactDOM_istep);
        reactDOM(:,istep)  = reactDOM_istep(:) ;
    end
    
    
    %%% STRESSES
    % ---------------
    DATAIN = DefaultField(DATAIN,'ISTWIST_ANGLE',0) ;
    if DATAIN.ISTWIST_ANGLE ==0
        T = RotationMatrix_z(ROTATION) ;
    else
        T = RotationMatrix_x(ROTATION) ;
    end
    
%     stressDOM = reshape(stressDOM,nstrain,[]) ;
%     stressDOM = T\stressDOM;
%     stressDOM  = stressDOM(:) ;
%     
       for istep = 1:size(stressDOM,2)
        stressDOM_istep = reshape(stressDOM(:,istep),nstrain,[]) ;
        stressDOM_istep = (T\stressDOM_istep);
        stressDOM(:,istep)  = stressDOM_istep(:) ;
    end
    
    
end





