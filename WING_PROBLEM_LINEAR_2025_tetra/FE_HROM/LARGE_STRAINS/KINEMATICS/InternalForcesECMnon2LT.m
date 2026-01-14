function [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2LT(OPERFE,PoneST,DATA,VAR,tauNONallDER_q) 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/InternalForcesECMnon2.m

nF = DATA.MESH.ndim^2 ;
% PoneST is the 1st Piola-Kirkhoff stresses for all the ECM points (master ECM points)
% PoneST = [PoneST_1; PoneST_2 ; \vdots PoneST_mECM_master]
% What we need is the internal work density matrix   \fINTredNON_{mst}
% such that
% (\FintREDnon{}{})_{mst} =   \fINTredNON_{mst} w_{mst}
% Here w_{mst} is OPERFE.wSTs
% The size of \fINTredNON_{mst} is:  as many columns as master ECM
% points
% What about the rows ?  Number of  generalized coordinates, including the constrained ones.


nmodesUall = size(OPERFE.Bst,2) ;
mMST =length(OPERFE.wSTs);  % Number of master integration points
fINTred_mst = zeros(nmodesUall,mMST);
for icomp = 1:nF
    icol = icomp:nF:size(PoneST,1) ;
    
    BstLOC =OPERFE.Bst(icol,:) ;
    % The size of the above matrix is
    %  mMST x nmodesUall
    PoneSTloc = PoneST(icol) ;
    % The size of the above matrix mMSTx1
    % Thus, we have to compute
    fINTloc = bsxfun(@times,BstLOC,PoneSTloc);
    % fINTloc is now  a mMST x nmodesUall matrix
    
    % Now we just add the contribution
    fINTred_mst = fINTred_mst + fINTloc';
end

% \item Project the above into the tangent space of the manifold:
%\fINTredNON_{mst,\DOFl} =  \tauNONder^T  \fINTred_{mst,\DOFl}
fINTredNON_mst = tauNONallDER_q'*fINTred_mst ;


% With   fINTredNON_mst at hand, we can compute the first contribution
% to the internal forces (that of the master point)
Fint  = fINTredNON_mst*OPERFE.wSTs ;
% NExt step: we use just the unconstrained component of FintREDnon_mst
fINTredNON_mst_DOFl = fINTredNON_mst(OPERFE.DOFl,:) ;
% as an input for the nonlinear mapping called "etanon"
% \fINTredNON_{slv,\DOFl} =\etaNON(\fINTredNON_{mst,\DOFl})
%    fINTredNON_slv_DOFl = OPERFE.etaNON(fINTredNON_mst_DOFl)' ;

[fINTredNON_slv_DOFl,der_fINTredNON_slv_DOFl,~] =...
    feval(OPERFE.DATA_regress_eta_der.nameFunctionEvaluate,fINTredNON_mst_DOFl,OPERFE.DATA_regress_eta_der) ;


% And then multiplied by their weights to get the contribution
% to the reduced internal forces
FintREDnon_slv_DOFl = fINTredNON_slv_DOFl*OPERFE.wRED_slv ;
% Now we add the slave contribution to the one from the master
% integration points
Fint(OPERFE.DOFl) = Fint(OPERFE.DOFl) + FintREDnon_slv_DOFl ;


% For linearization purposes, we need to compute the  the ``equivalent'' integration weight
% for the master point
%
%  \wwMSTeq \defeq   w_{mst} + \etaNONder(\fINTredNON_{mst,\DOFl}) \w_{mst}

%wMSTeq = OPERFE.wSTs + OPERFE.etaNONder(fINTredNON_mst_DOFl)'*OPERFE.wRED_slv ;

wMSTeq = OPERFE.wSTs + der_fINTredNON_slv_DOFl'*OPERFE.wRED_slv ;




%  Compute the extended internal forces at the  master point assuming the above is
%  the integration weight
%   \FintREDeqMSTl \defeq  \fINTred_{mst,\DOFl}    \wwMSTeq
FintREDeqMST = fINTred_mst*wMSTeq ;
% We use this to define a residual (that will be employed for constructing
% the reduced stiffness matrix)
%  Calculate the corresponding residual (Eq. \ref{eq:reis33}):
% \ResREDeqMSTl \defeq \FintREDeqMSTl-\FextRED{}{\DOFl}
ResREDeqMST = FintREDeqMST-VAR.FEXT;
