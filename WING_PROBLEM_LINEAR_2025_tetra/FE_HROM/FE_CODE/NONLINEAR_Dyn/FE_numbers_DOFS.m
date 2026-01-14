function  OPERfe = FE_numbers_DOFS(nnode,ndim,DOFs,DOFm,Gb,nelem,DATA,nnodeE,ngausT)


% Set of nodes at which displaceemnts are unknown
DATA= DefaultField(DATA,'ndof',nnode*ndim) ;  
DOFl = (1:DATA.ndof)' ;  % All DOFs
DOFl([DOFs ;DOFm]) = [] ; % All DOFs but slave DOFs (or DOFr) and master (DOFm)
DOFf = [DOFl ; DOFm] ;   % All unknowns, first interior DOFs, and then master DOFs 
% --------------------------
OPERfe.DOFf = DOFf ; % Change of notation !!!!  ---> See BeamROM_nonlinear
OPERfe.DOFl = DOFl ; 
OPERfe.DOFm = DOFm ; 
OPERfe.DOFs = DOFs ; 
OPERfe.Gbound = Gb ;
OPERfe.nnode = nnode ;
OPERfe.ndim = ndim ;
OPERfe.nelem = nelem ;
OPERfe.nnodeE = nnodeE ;
OPERfe.nstrain = DATA.nstrain ; 
OPERfe.ngaus = ngausT/OPERfe.nstrain ;  % Total number of Gauss points
OPERfe.ngausE = OPERfe.ngaus/nelem ; 

% ------------------
OPERfe.ndof = DATA.ndof; 
OPERfe.ngausT = ngausT  ;%*nstrain ;
% ------------------------------------------
OPERfe.nsteps = length(DATA.TIME_DISCRETIZATION);  % Number of time steps
