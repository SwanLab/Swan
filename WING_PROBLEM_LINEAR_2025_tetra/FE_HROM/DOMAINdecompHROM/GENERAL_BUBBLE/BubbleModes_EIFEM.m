function GammaBUB = BubbleModes_EIFEM(PsiSEf,PhiDEF,BasisUdeform,b,PhiRB,Mdom,DATALOC,MESH,MdomCHOL,DATAoffline,MintfINV_chol)
%==========================================================================
% function GammaBUB = BubbleModes_EIFEM(...)
%
% Computes the bubble modes ΓBUB used in the EIFEM (Empirical Interscale
% Finite Element Method) framework, based on a hybrid approach combining
% kinetic and zero-work definitions of bubble modes.
%
% AUTHOR:
% Joaquín Alberto Hernández Ortega
% Balmes 185, Barcelona
%
% DATE:
% Modified implementation on 18-05-2025 (Sunday)
%
% REFERENCES:
% - Appendix "How to compute bubble modes" in EIFEM_largeROTfinal.pdf
% - BubbleModes_BasicSEfromDEF.m and DeformModesNONL_INVAR.m
%
% THEORY SUMMARY:
% The method computes a set of bubble modes Γ such that:
%   (i) span([ΦDEF, Γ]) = span(Υ), with Υ = [ΦDEF, Φin]
%   (ii) Γb ≈ orthogonal to ΨSEf (interface self-equilibrated modes)
%
% Two types of bubble modes are computed:
% - Kinematic bubble modes (Γ_kin): modes with vanishing interface trace,
%   obtained by projecting Υ onto the complement of D = Υ(V_d * S_d^{-1}).
% - Zero-work bubble modes (Γ_w): when p = dim(D) - dim(ΦDEF) > 0, the
%   subspace of D most orthogonal to ΨSEf is extracted via SVD and used
%   to define Γ_w with minimal interaction (work) with ΨSEf.
%
% The final matrix ΓBUB combines both Γ_kin and Γ_w, and is M-orthogonalized.
%
% OPTIONS & NOTES:
% - Purging of rigid body modes is available but disabled by default.
% - Compatible with both direct and Cholesky-based SVDs.
% - Principal angles are displayed to check the linear independence between
%   ΓBUB and ΦDEF, and the degree of orthogonality with ΨSEf.
%
%==========================================================================


if nargin == 0
    load('tmp1.mat')
end

disp(['Method based on orthogonal projections  '])
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/06_ForcingBubbles.mlx
%  /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf

%\item  Let us begin by    introducing an alternative decomposition of the deformational space;
%  these decomposition involves the strain modes $\PhiDEF$ and a set of modes $\GammaBUB$
% (linearly independent from $\PhiDEF$) such that, on the one hand
%\begin{equation}
% \spanb{[\PhiDEF, \GammaBUB]}  = \spanb{\UpsilonDEF}
%\end{equation}
%where   $\UpsilonDEF = [\PhiDEF,\PhiDEFin]$, and  on the other,
% $\GammaBUBb =\GammaBUB(\bDOFS,:) $ is ``as orthogonal'' as possible to $\PsiSE$.
UpsilonDEF = [PhiDEF,BasisUdeform.PhiDEFcomp] ; 
% The latter condition naturally raises the question on how to mathematically state the condition
%``as orthogonal as possible to $\PsiSE$''.     Consider the SVD of the interface entries of the deformational modes $\UpsilonDEF$:
%\begin{equation}
%[\Dbar,\S_d,\V_d] \leftarrow \SVD{\UpsilonDEFb,\epsilon_d}
%\end{equation}
%where tolerance $\epsilon_d$ is the truncation tolerance
%  (it   may be set to $\epsilon_d \approx 10^{-6}$).
epsilon_d = 0;
DATALOC.RELATIVE_SVD = 1;
[Dbar,S_d,V_d] = SVDT(UpsilonDEF(b,:),epsilon_d,DATALOC) ;
% It follows from the definition of SVD that, for the boundary DOFs:
%     \begin{equation}
%   \Dbar = \UpsilonDEFb (\V_p \S^{-1})
%  \end{equation}
%  and, therefore,  for all the DOFs
%   \begin{equation}
%   \D = \UpsilonDEF (\V_p \S^{-1})
%  \end{equation}
VSd = bsxfun(@times,V_d',1./S_d)' ;
D = UpsilonDEF*VSd ;
GammaBUB_kin = [] ;
if size(D,2) ~= size(UpsilonDEF,2)
    % New snippet (implemented on May 18th 2025), see problem in
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN/ONE_CELL_COMPRESSION.mlx
    % \textbf{Genuine bubble modes (zero at the interface DOFs)}.  If $\ncol{\D} < \ncol{\UpsilonDEFb}$,
    % then it is inferred that the first $\ncol{\UpsilonDEFb} - \ncol{\D}$ left singular vectors of
    %\begin{equation}
    % \J =  \UpsilonDEF - \D (\D^T \D)^{-1} \D^T \UpsilonDEF
    % \end{equation}
    coeff = (D'*Mdom*D)\(D'*Mdom*UpsilonDEF) ;
    J = UpsilonDEF-D*coeff ;
    nkinetic =  size(UpsilonDEF,2)-size(D,2)   ;
    if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
        %  if ~isempty(MdomCHOL)
        DATALOCww.Mchol = MdomCHOL ;
        [GammaBUB_kin,~,~] = WSVDT(J,[],DATALOCww) ;
        % else
        %     [GammaBUB_kin,~,~] = WSVDT(J,Mdom) ;
        % end
    else
        [GammaBUB_kin,~,~] = SVDT(J) ;
    end
    GammaBUB_kin = GammaBUB_kin(:,1:nkinetic) ;
    
end
p  = size(D,2) - size(PhiDEF,2) ;
[P,SS1,VV] = SVDT(PsiSEf) ;
% Projection of Dbar onto the orthogonal complement of span(P)
T = Dbar- P*(P'*Dbar);
[G,St,Vt] = SVDT(T) ;
[Jbar, SJ,VJ] = SVDT(G'*Dbar) ;
disp(['Principal Angles formed by the kernel of PsiSE and the columnspace of UpsilonDEF'])
disp(['If the first ',num2str(p),'angle(s) are zero, it means the bubble condition is exactly met (no need for approximations)'])
real(acosd(SJ))
Y = D*VJ(:,1:p) ;

s = 1:size(Y,1) ;
s(b)  = [] ;
GammaBUB = zeros(size(Y,1),p) ;
GammaBUB(s,:) = Y(s,:) ;
GammaBUB(b,:) = G*Jbar(:,1:p);


GammaBUB = [GammaBUB_kin,GammaBUB] ;

if  DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
    %  [GammaBUB,SS,VV] = WSVDT(GammaBUB,Mdom) ;
    DATALOCww.Mchol = MdomCHOL ;
    [GammaBUB,~,~] = WSVDT(GammaBUB,[],DATALOCww) ;
    
else
    [GammaBUB,SS,VV] = SVDT(GammaBUB,[]) ;
end


disp('Checking that GammaBUB and PhiDEFf are linearly independent')
disp('If this is so, all the angles formed by the corresponding subspaces should greater than 0') ;
[uuu,sss,vvv] = SVDT(GammaBUB'*Mdom*PhiDEF) ;
alpha = real(acosd(sss))
% Let
if any(alpha < 1.0)  %%%The
    warning('Some bubble modes are almost linearly dependent to the basic modes')
end
disp(['Ckecking that the modes are indeed free-residual bubble modes '])


PURGING_RIGID_BODY_MODES = 0;
if  PURGING_RIGID_BODY_MODES == 1
    % Option introduced on April 24th 2025, see
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/10_AUXETIC_2D.mlx
    % This does not alter the property GammaBUB(b,:)'*PsiSEf = zero, because  (GammaBUB(b,:)' + A*R')*PsiSEf = zero
    
    % However, strictly speaking, in large strains, only TRANSLATIONS, are
    % to be purged
    
    ONLY_translations = 1;
    if ONLY_translations == 1
        % This should be the default option for large rotations
        if size(PhiRB,2)  == 3
            nTR = 1:2 ;
        elseif size(PhiRB,2)  == 6
            nTR =1:3 ;
        end
        GammaBUB = SprojDEF_operator(PhiRB(:,nTR),Mdom,GammaBUB) ;
        
    else
        
        GammaBUB = SprojDEF_operator(PhiRB,Mdom,GammaBUB) ;
    end
    
end


if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
    
    DATALOCww.Mchol = MintfINV_chol ;
    [Z,SS,VV] =WSVDT(GammaBUB(b,:),[],DATALOCww) ;  % Z*SS*VV^T = Ub  --> U(b,:)*VV*inv(SS) = Z
else
    [Z,SS,VV] =WSVDT(GammaBUB(b,:),[]) ;
    
end


[uu,ss,vv] = SVDT(Z'*PsiSEf) ;
disp('-------------------------------------------------')
disp('PRINCIPAL  ANGLES between  PsiSE     AND  GammaBUB(b,:) )')
disp('-------------------------------------------------')
theta = real(acosd(ss))
DATALOC.NAME_BASE = 'BUBBLE';
PlotModesDEF_SubdomainLevel(DATALOC,GammaBUB,MESH);

