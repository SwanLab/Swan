function [stressDOM,reactDOM,Ndom,wDIAG_RHS,Bdom,Wdom,K,CgloDOM] = ...
    ReactionsDOMnonl(Nst,fNOD,ListGauss,ndim,DOFs,Fpnt,CNb,NODESrve,Tnod,...
    TypeElementB,CONNECTb,...
    Bst,wSTs_RHS,wSTs,COOR,nnodeDOM,faceDOFS,BasisUrb,nstrain,idom,iproject,...
    DATAIN,stressGLO,DATA_INPUT_FE,Cglo)
% Computing reactions and stresses at each domain 
if nargin == 0
    load('tmp.mat')
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

% stresses 
stressDOM = stressGLO(ListGaussDOFS,:) ; 
% Internal forces 
fintDOM = (wDIAG*Bdom)'*stressDOM ;

nsteps = size(fintDOM,2) ; 

% Elasticity matrix 
CgloDOM = [] ; 
if  ~isempty(Cglo)
CgloDOM=  Cglo(ListGaussDOFS,ListGaussDOFS) ; 
else
     load(DATAIN.NAME_WS_MODES,'CgloDOM') ;
end

 K = [] ;  % Stiffness matrix
if idom == 1 && iproject==1
    K = (Cglo(ListGaussDOFS,ListGaussDOFS)*Bdom) ;
    K = (wDIAG*Bdom)'*K ;
end

%%%% Reactions
% See function ExtForces_Pdisp_TIME.m
ftracDOM = bsxfun(@times,ftracDOM',DATA_INPUT_FE.FACTOR_TIME_TRACTION_FORCES')' ;
ftracDOM = ftracDOM(:,1:nsteps) ;

if length(fbodyDOM) >1
    fbodyDOM = bsxfun(@times,fbodyDOM',DATA_INPUT_FE.FACTOR_TIME_BODY_FORCES')' ;
    fbodyDOM = fbodyDOM(:,1:nsteps) ;
end


reactDOMorig = fintDOM - (ftracDOM+ fbodyDOM ) ;


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
reactDOM(faceDOFS,:) = reactDOMorig(faceDOFS,:)  - BasisRrb(faceDOFS,:)*(BasisRrb(faceDOFS,:)\reactDOMorig(faceDOFS,:)) ;
