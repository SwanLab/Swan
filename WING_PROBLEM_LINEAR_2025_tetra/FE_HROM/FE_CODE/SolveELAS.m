function [d strainGLO stressGLO  React posgp DATAOUT] = ...
    SolveELAS(K,Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,celasglo,...
    typePROBLEM,celasgloINV,Cglo,Bst,DOFm,Gb,DATA,wST,Nst) ;
% This function returns   the (nnode*ndim x 1) vector of nodal displacements (d),
% as well as the arrays of stresses and strains
%%% points  (qheatGLO)
% Input data
% K = Global stiffness matrix   (nnode*ndim x nnode*ndim)
% Fb = External force vector due to  body forces  (nnode*ndim x 1)
% Ftrac = External force vector due to  boundary tractions    (nnode*ndim x 1)
% DOFr = Set of restricted DOFs
% dR = Vector of prescribed displacements
% ----------------------
%dbstop('14')
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;     %
DATA.nstrain = size(celasglo,1) ;
% Solution of the system of FE equation
% Right-hand side
F = Fb + Ftrac ;
% Set of nodes at which temperature is unknown
DOFf = (1:nnode*ndim)' ;
DOFf([DOFr ;DOFm]) = [] ;
DOFl = [DOFm ; DOFf] ;




% dL =  K^{-1}*(Fl .Klr*dR)
%dbstop('28')

if isempty(Gb)
    Kll =  K(DOFl,DOFl) ;
    Fl = F(DOFl) ;
    disp(['length(DOFl) =',num2str(length(DOFl))]) ;
    disp(['length(DOFr) =',num2str(length(DOFr))]) ;
    disp(['length(dR) =',num2str(length(dR))]) ;
    KlrU = K(DOFl,DOFr)*dR ;
else
    
    
    Gast = [Gb,sparse(size(Gb,1),length(DOFf))] ;
    Kll = K(DOFl,DOFl) + (Gast'*K(DOFr,DOFl) + K(DOFl,DOFr)*Gast) + Gast'*(K(DOFr,DOFr)*Gast);
    Fl = [Gb'*F(DOFr) + F(DOFm); F(DOFf)] ;
    KlrU = (Gast'*K(DOFr,DOFr)+K(DOFl,DOFr))*dR ;
    % Klr = (Gast'*K(DOFr,DOFr)+K(DOFl,DOFr)) ;  % BORRAR
    
end



%dbstop('48')
disp('------------------------------------')
disp('SOLVING  dL = inv(Kll)* FL')
disp('------------------------------------')
aa = tic;
if DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 1 || DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 3
    FreeingMemory ;  % Deleting Bst, Cglo....
elseif DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 2
    FreeingMemoryAndStoreSomeOfThem ;  % Deleting Bst, Cglo....
end

BBBB = whos('Kll') ;
disp(['Size Kll = ',num2str(BBBB(1).bytes*1e-9),' Gb' ])
if DATA.TYPESOLVER  == 0
    %     pp = symrcm(Kll) ;
    %     spy(Kll(pp,pp));
    disp('Direct solver ...')
    
    dL =Kll\(Fl-KlrU) ;
    
elseif DATA.TYPESOLVER  >=1
    dL =  IterativeSolverFE(DATA,Kll,Fl,KlrU) ; % Iteratives solvers
end
%dL =Kll\(Fl-KlrU) ;
aa = toc(aa) ;
disp(['DONE in:  ',num2str(aa),' s']);
disp('------------------------------------')
disp('------------------------------------')

% Vector of   displacements
%dbstop('44')
d = zeros(nnode*ndim,1) ;
if isempty(Gb)
    d(DOFl)= dL ;
    d(DOFr) = dR ;
else
    d(DOFl)= dL ;
    d(DOFr) = dR + Gb*dL(1:length(DOFm));  %dl
end
% Reaction forces
if DATA.PLOT.REACTIONS == 1
    React = zeros(size(d)) ;
    React(DOFr) = K(DOFr,:)*d -F(DOFr) ;
    if ~isempty(Gb)
        React(DOFm) = - Gb'*React(DOFr) ;
    end
else
    React = [] ;
end

%%%% COmputation of heat flux vector at each gauss point
disp('Computation of stress and strains at each Gauss point')

if DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING ==  1
    disp('Stresses will not be printed in GID (because of option DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING =  1)')
    strainGLO =[]; stressGLO =[] ; posgp = []; DATAOUT = [] ;
else
    
    if DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 2
        disp('loading Bst ....')
        load(NWS_Bst),
        disp('loading Cglo')
        load(NWS_Cglo) ;
        celasglo = []; ;celasgloINV=[] ; Nst = [] ;
    end
    
    [strainGLO stressGLO posgp DATAOUT ]= ...
        StressStrains(COOR,CN,TypeElement,celasglo,d,typePROBLEM,celasgloINV,...
        Bst,Cglo,DATA,wST,Nst,React) ;
    
    if  DATA.MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING == 2
        delete(NWS_Bst)
        delete(NWS_Cglo)
    end
    
end


DATAOUT.BOUNDDATA.DOFm = DOFm ;
DATAOUT.BOUNDDATA.DOFr = DOFr ;
DATAOUT.BOUNDDATA.DOFf = DOFf ;
DATAOUT.BOUNDDATA.Gb   = Gb ;
DATAOUT.Kll   = [] ;  % BORRAR
DATAOUT.Klr   = [] ;  % BORRAR
save(DATA.nameWORKSPACE,'DOFr','-append')
