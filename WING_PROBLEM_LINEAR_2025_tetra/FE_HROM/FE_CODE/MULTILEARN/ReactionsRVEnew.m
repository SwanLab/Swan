function [reactDOM stressDOM,Bdom,Wdom] = ReactionsRVEnew(DATA,NODESfaces,dRVE,...
    nnode,ndim,Nst,wSTs_RHS,fNOD,NODESrve,CNrve,posgp,IndElements,posgp_RHS,...
    CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,wSTs,Bst,TypeElement,BasisRrb,Cglo,d,...
    DATACUBATURE,DATAIN)
% See DomainDecom_SVD.m
% Computing self-equilibrated reaction forces
% --------------------------------------------
% -----------------------
%dbstop('10')
if nargin == 0
    load('tmp.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing body forces at domain idom
% -------------------------------------
%dbstop('16')
fbodyDOM = 0 ;
if ~isempty(fNOD) & any(fNOD)
    % Set of Gauss points associated to domain idom
    % -----------------------------------------------------------------------
    ngaus = size(posgp_RHS,2) ; % Number of Gauss points
    nelem = length(IndElements) ;
    ListGauss =  small2large(IndElements,ngaus) ;  % List of GAuss points
    % Set of DOFs associated to domain idom
    % ---------------------------------------------
    ndim = size(COORabs,2) ; nnode = length(NODESrve) ;
    ListDOFs =  small2large(NODESrve,ndim) ;  % List of DOFs
    % Therefore
    ListGaussDOF =  small2large(ListGauss,ndim) ;
    Ndom = Nst(ListGaussDOF,ListDOFs) ;
    WdomN = wSTs_RHS(ListGauss)  ;
    fbodyNODES = fNOD(ListDOFs) ;
    % Converging WdomN into a diagonal matrix affecting all DOFs
    wDIAG = CompWeightDiag(WdomN,ndim)  ;
    fbodyDOM = (wDIAG*Ndom)'*(Ndom*fbodyNODES) ;
else
    ngaus = size(posgp_RHS,2) ; % Number of Gauss points
    nelem = length(IndElements) ;
    ListGauss =  small2large(IndElements,ngaus) ;  % List of GAuss points
    ListDOFs =  small2large(NODESrve,ndim) ;  % List of DOFs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing traction forces -->  FtracDOM
% ------------------------------------------------------
% WE shall invoke function
%  Ftrac = FtracCOMPvect(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb) ;
% but adapting the inputs for domain idom
FpntDOM = zeros(size(Fpnt)) ;
FpntDOM(ListDOFs) = Fpnt(ListDOFs)  ; % Point loads

CNbDOM = cell(size(CNb)) ;
TnodDOM = cell(size(CNbDOM)) ;
%dbstop('57')
for idim = 1:ndim
    if ~isempty(CNb{idim})
        [CNbDOM{idim} ListEbndDOM]= ElemBnd(CNb{idim},NODESrve) ;
        TnodDOM{idim} = Tnod{idim}(ListEbndDOM,:) ;
    end
end
ftracDOM = FtracCOMPvect(COOR,CNbDOM,TypeElementB,FpntDOM,TnodDOM,CONNECTb) ;
ftracDOM = ftracDOM(ListDOFs) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of internal forces
Wdom =  wSTs(ListGauss);
if ndim == 2
    nstrain = 3 ;
elseif ndim ==3
    nstrain = 6 ;
end
%dbstop('77')
wDIAG = CompWeightDiag(Wdom,nstrain)  ;
ListGaussDOFS =  small2large(ListGauss,nstrain) ;
Bdom = Bst(ListGaussDOFS,ListDOFs) ;
%load('tmp0.mat')
%dbstop('82')
% Computing stresses and internal forces
% ------------------------
strainDOM = Bdom*dRVE ;
%strainDOM = Bdom*d(ListDOFs) ;

stressDOM = (Cglo(ListGaussDOFS,ListGaussDOFS)*strainDOM) ; %
for istrain =1:nstrain
    stressDOM(istrain:nstrain:end) = stressDOM(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
end
%dbstop('90')
fintDOM = (wDIAG*Bdom)'*stressDOM ;



%%%% Reactions
reactDOMorig = fintDOM - (ftracDOM+ fbodyDOM) ;
%%% Reaction vectors at faces f1 and f2
% f1 = > NODESfaces{idom}{3}
% f2 = > NODESfaces{idom}{1}
%dbstop('100')
%DATAIN.TYPE_REPETITION = 'MULTIDIRECTION';
DATAIN = DefaultField(DATAIN,'REPETITION_SEVERAL_DIRECTION',0) ; 
if  DATAIN.REPETITION_SEVERAL_DIRECTION == 0
    %-------------------
    % 1-D tiled problems
    % ------------------
    f1 = NODESfaces{1} ;
    f2 = NODESfaces{3} ;
    f1 = small2large(f1,ndim) ;
    f2 = small2large(f2,ndim) ;
    f = [f1; f2] ;
else
    % In 2D-tiled and 3D-tiled problems, we shall include the whole set of
    % Boundary nodes in set "f"
  %  dbstop('112')
    f = [] ;
    for ifaces = 1:length(NODESfaces)
        f = [f; NODESfaces{ifaces}'] ;
    end
    f = unique(f) ; 
    f = small2large(f,ndim) ;
end


% SELF-EQUILIBRATED MODES
% -----------------------
% %%% We seek self-equilibrated reaction force vectors
% Interior nodes are   not accounted for9
DOFtot = 1:ndim*nnode ;
DOFi = setdiff(DOFtot,f) ;
BasisRrb(DOFi,:) = 0 ;
% Reaction forces with resultant equal to zero
% ---------------------------------------------
%dbstop('145')
reactDOM = zeros(size(reactDOMorig)) ;
reactDOM(f) = reactDOMorig(f)  - BasisRrb(f,:)*(BasisRrb(f,:)\reactDOMorig(f)) ;




%
%
% PRUEBAS_PRINT_MODES = 0;
% %dbstop('138')
% if PRUEBAS_PRINT_MODES ==1
%     DATA.NODES = NODESrve' ;
%
%     GidPostProcess(COOR,CNrve{idom},TypeElement,[],[], ...
%         [],  reactDOMorig,[''],posgp,['prueba.msh'],[],DATA);
% end



