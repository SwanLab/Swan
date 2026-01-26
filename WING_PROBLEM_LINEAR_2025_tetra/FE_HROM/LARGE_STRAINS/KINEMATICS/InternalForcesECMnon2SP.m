function [Fint,FintNON,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] =...
    InternalForcesECMnon2SP(OPERFE,PK2STRESS_incre,DATA,VAR,tauNONallDER_q) 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/InternalForcesECMnon2.m
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECM.mlx

if nargin == 0
    load('tmp2.mat')
end

nF = DATA.MESH.nstrain  ;
% PK2STRESS_incre is the  nonlinear stresses for all the ECM points (master ECM points)
 % What we need is the nonlinear component of internal work density matrix   \fINTredNON_{mst}
% such that
% (\FintNONREDnon{}{})_{mst} =   \fINTredNON_{mst} w_{mst}
% Here w_{mst} is OPERFE.wSTs
% The size of \fINTredNON_{mst} is:  as many columns as master ECM
% points
% What about the rows ?  Number of  generalized coordinates, including the constrained ones.


nmodesUall = size(OPERFE.Bst,2) ;
mMST =length(OPERFE.wSTs);  % Number of master integration points
fINTredNON_mst = zeros(nmodesUall,mMST);
for icomp = 1:nF
    icol = icomp:nF:size(PK2STRESS_incre,1) ;
    
    BstLOC =OPERFE.Bst(icol,:) ;
    % The size of the above matrix is
    %  mMST x nmodesUall
    PK2STRESS_increloc = PK2STRESS_incre(icol) ;
    % The size of the above matrix mMSTx1
    % Thus, we have to compute
    fINTloc = bsxfun(@times,BstLOC,PK2STRESS_increloc);
    % fINTloc is now  a mMST x nmodesUall matrix
    
    % Now we just add the contribution
    fINTredNON_mst = fINTredNON_mst + fINTloc';
end

% \item Project the above into the tangent space of the manifold:
%\fINTredNON_{mst,\DOFl} =  \tauNONder^T  \fINTred_{mst,\DOFl}
fINTredNON_mst = tauNONallDER_q'*fINTredNON_mst ;


% With   fINTredNON_mst at hand, we can compute the first contribution
% to the internal forces (that of the master point)
FintNON  = fINTredNON_mst*OPERFE.wSTs ;

% NExt step: we use just the unconstrained component of FintNONREDnon_mst
fINTredNON_mst_DOFl = fINTredNON_mst(OPERFE.DOFl,:) ;
% as an input for the nonlinear mapping called "etanon"
 

fINTredNON_slv_DOFl = zeros(size(fINTredNON_mst_DOFl,1),length(OPERFE.wRED_slv)) ; 
% We loop over the number of rows   of fINTredNON_slv_DOFl
for icomp = 1:length(OPERFE.DATA_regress_eta_der)
[fINTredNON_slv_DOFl(icomp,:),der_fINTredNON_slv_DOFl,~] =...
    feval(OPERFE.DATA_regress_eta_der{icomp}.nameFunctionEvaluate,fINTredNON_mst_DOFl(icomp),OPERFE.DATA_regress_eta_der{icomp}) ;
end
% And then multiplied by their weights to get the contribution
% to the reduced internal forces
FintNONREDnon_slv_DOFl = fINTredNON_slv_DOFl*OPERFE.wRED_slv ;
% Now we add the slave contribution to the one from the master
% integration points
FintNON(OPERFE.DOFl) = FintNON(OPERFE.DOFl) + FintNONREDnon_slv_DOFl ;

% Now we add the elastic part 
Fint =  tauNONallDER_q'*(OPERFE.KstiffLINEAR*VAR.DISP) +  FintNON ;  


% For linearization purposes, we need to compute the  the ``equivalent'' integration weight
% for the master point
%
%  \wwMSTeq \defeq   w_{mst} + \etaNONder(\fINTredNON_{mst,\DOFl}) \w_{mst}

%wMSTeq = OPERFE.wSTs + OPERFE.etaNONder(fINTredNON_mst_DOFl)'*OPERFE.wRED_slv ;

wMSTeq = OPERFE.wSTs + der_fINTredNON_slv_DOFl'*OPERFE.wRED_slv ;




%  Compute the extended internal forces at the  master point assuming the above is
%  the integration weight
%   \FintNONREDeqMSTl \defeq  \fINTred_{mst,\DOFl}    \wwMSTeq
FintNONREDeqMST = fINTredNON_mst*wMSTeq ;
% We use this to define a residual (that will be employed for constructing
% the reduced stiffness matrix)
%  Calculate the corresponding residual (Eq. \ref{eq:reis33}):
% \ResREDeqMSTl \defeq \FintNONREDeqMSTl-\FextRED{}{\DOFl}
ResREDeqMST = FintNONREDeqMST-VAR.FEXT;
