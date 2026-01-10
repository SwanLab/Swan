function [reactDOM] = ReactionsRVE_conflict(idom,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
    COORrve,Nst,wSTs_RHS,fNOD,NODESrve,CNrve,posgp,IndElements,posgp_RHS,...
    CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,stressGLO,wSTs,Bst,TypeElement)
% See DomainDecom_SVD.m
% Computing self-equilibrated reaction forces
% --------------------------------------------
% -----------------------
if nargin == 0
    load('tmp1.mat')
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
    ndim = size(COORabs,2) ; nnode = length(NODESrve{idom}) ;
    ListDOFs =  small2large(NODESrve{idom},ndim) ;  % List of DOFs
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
    
    ndim = size(COORabs,2) ; nnode = length(NODESrve{idom}) ;
    
    ListDOFs =  small2large(NODESrve{idom},ndim) ;  % List of DOFs
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
        [CNbDOM{idim} ListEbndDOM]= ElemBnd(CNb{idim},NODESrve{idom}) ;
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
wDIAG = CompWeightDiag(Wdom,nstrain)  ;
ListGaussDOFS =  small2large(ListGauss,nstrain) ;
Bdom = Bst(ListGaussDOFS,ListDOFs) ;
stressDOM = stressGLO(ListGaussDOFS) ;
%dbstop('79')
fintDOM = (wDIAG*Bdom)'*stressDOM ;

%%%% Reactions
reactDOMorig = fintDOM - (ftracDOM+ fbodyDOM) ;
%%% Reaction vectors at faces f1 and f2
% f1 = > NODESfaces{idom}{3}
% f2 = > NODESfaces{idom}{1}
f1 = NODESfaces{idom}{1} ;
f2 = NODESfaces{idom}{3} ;
f1 = small2large(f1,ndim) ;
f2 = small2large(f2,ndim) ;


f = [f1; f2] ;


% SELF-EQUILIBRATED MODES
% -----------------------
%%% We seek self-equilibrated reaction force vectors

%   Therefore, we have to construct the rigid body modes ()
COORrel = COORrve{idom};  % Coordinates with respect to the reference point
if ndim == 2
    BasisRrb =  zeros(prod(size(COORrel)),3) ;
    BasisRrb(1:ndim:end,1) = 1;
    BasisRrb(2:ndim:end,2) = 1;
    BasisRrb(1:ndim:end,3) = -COORrel(:,2);
    BasisRrb(2:ndim:end,3) = COORrel(:,1);
else
    
    BasisRrb =  zeros(prod(size(COORrel)),6) ;
    BasisRrb(1:3:end,1) = 1 ;
    BasisRrb(2:3:end,2) = 1;
    BasisRrb(3:3:end,3) = 1;
    % Rotation modes
    BasisRrb(2:3:end,4) = COORrel(:,3);
    BasisRrb(3:3:end,4) = -COORrel(:,2);
    
    BasisRrb(1:3:end,5) = -COORrel(:,3);
    BasisRrb(3:3:end,5) = COORrel(:,1);
    
    BasisRrb(1:3:end,6) = COORrel(:,2);
    BasisRrb(2:3:end,6) = -COORrel(:,1);
    
end

% Interior nodes are   not accounted for
DOFtot = 1:prod(size(COORrel)) ;
DOFi = setdiff(DOFtot,f) ;
BasisRrb(DOFi,:) = 0 ;
% Reaction forces with resultant equal to zero
% ---------------------------------------------
%dbstop('132')
reactDOM = zeros(size(reactDOMorig)) ;
reactDOM(f) = reactDOMorig(f) - BasisRrb(f,:)*(BasisRrb(f,:)\reactDOMorig(f)) ;






PRUEBAS_PRINT_MODES = 0;
%dbstop('138')
if PRUEBAS_PRINT_MODES ==1
    DATA.NODES = NODESrve{idom}' ;
    
    GidPostProcess(COOR,CNrve{idom},TypeElement,[],[], ...
        [],  reactDOMorig,[''],posgp,['prueba.msh'],[],DATA);
end



