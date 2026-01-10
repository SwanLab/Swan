function [mNORMst,kM] = NormalNonUnitaryBoundary(OPERFEH,VAR,DATA)
% See formulation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
% Normal (nonunitary) vector at each Gauss point (as well as its contribution to the tangent matrix )
%
%  \begin{equation}
%  \boxed{\mNORMst = \mNORMiniST  + \Par{\dmNORMxdST  +  \dmNORMdxST  + \dmNORMddST}}
% \end{equation}

% % wHERE
% %  \begin{equation}
% \label{eq:8ddg-}
%  \dmNORMxdST^e  =  \tTANGiniST{e}_1 \times  \BstB{2}{e}   \d^e
% \end{equation}

% \begin{equation}
%  \label{eq:::8}
%  \kNORMxdST^e \defeq \derpar{\dmNORMxdST^e}{\d^e}  = [\tTANGiniST{e}_1 \times  \BstB{2}{e}]_{3 \ngausFT \times 3 \nnodeEb{e}{}}
% \end{equation}

% In summary **** m = mINI + k_xd*d  +  k_dx*d    +  k_dd(d)*d
% where
%  k_xd = tau_1_ini x  B_2
%  k_dx = -tau_2_ini x B_1
%  k_dd = (B_1*d) x B_2  -  (B_2*d) x B_1

% JAHO- 1-July-2021


ndim = 3;
ndimB = 2;

% -----------------------------------------------------------
% kM = dm/dd  * Contribution of m to the tangent matrix
% -----------------------------------------------------------
nrows = size(OPERFEH.mNORMiniST,1) ;
ncols = size(OPERFEH.BstB{1},2) ;

if  DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO ==1
    tau_1 = OPERFEH.tTANGiniST{1} + ConvertBlockDiag_general( OPERFEH.BstB{1},ndim,OPERFEH.irowsNDIM,OPERFEH.icolsNDIM)*OPERFEH.Lbool*VAR.DISP(:) ;
    tau_2= OPERFEH.tTANGiniST{2} +  ConvertBlockDiag_general( OPERFEH.BstB{2},ndim,OPERFEH.irowsNDIM,OPERFEH.icolsNDIM)*OPERFEH.Lbool*VAR.DISP(:) ;
    
    %     tau_1 = OPERFEH.tTANGiniST{1} + dtau_1 ;
    %         tau_2= OPERFEH.tTANGiniST{2} + dtau_2 ;
    %
    tau1 = reshape(tau_1,ndim,[]) ;
    tau2 = reshape(tau_2,ndim,[]) ;
    mNORMst =  cross(tau1,tau2) ;
    mNORMst = mNORMst(:);
    kM = [] ;
    %         EEE = norm(mNORMst-mNORMst_Alt)/norm(mNORMst_Alt) ;
    %         disp(['ERROR METHOD = ',num2str(EEE)])
    
else
    
    kM = zeros(nrows,ncols) ;
    %
    for icol = 1:ncols
        %   k_xd = tau_1_ini x  B_2
        k_xd_LOC = cross(reshape(OPERFEH.tTANGiniST{1},ndim,[]),reshape(OPERFEH.BstB{2}(:,icol),ndim,[])) ;
        %  k_dx = -tau_2_ini x B_1
        k_dx_LOC = -cross(reshape(OPERFEH.tTANGiniST{2},ndim,[]),reshape(OPERFEH.BstB{1}(:,icol),ndim,[])) ;
        
        % K_dd_1 = dtau_1 x B_2
        % dtau_1 = (B_1*d)
        dtau_1 = ConvertBlockDiag_general( OPERFEH.BstB{1},ndim,OPERFEH.irowsNDIM,OPERFEH.icolsNDIM)*OPERFEH.Lbool*VAR.DISP(:) ;
        k_dd_1_LOC = cross(reshape(dtau_1,ndim,[]),reshape(OPERFEH.BstB{2}(:,icol),ndim,[])) ;
        
        % K_dd_2 = -dtau_2 x B_1
        % dtau_2 = (B_2*d)
        dtau_2= ConvertBlockDiag_general( OPERFEH.BstB{2},ndim,OPERFEH.irowsNDIM,OPERFEH.icolsNDIM)*OPERFEH.Lbool*VAR.DISP(:) ;
        k_dd_2_LOC = - cross(reshape(dtau_2,ndim,[]),reshape(OPERFEH.BstB{1}(:,icol),ndim,[])) ;
        
        kLOC = k_xd_LOC + k_dx_LOC + k_dd_1_LOC   + k_dd_2_LOC ;
        
        kM(:,icol) = kLOC(:) ;
        
    end
    mNORMst = OPERFEH.mNORMiniST + ConvertBlockDiag_general( kM,ndim,OPERFEH.irowsNDIM,OPERFEH.icolsNDIM)*OPERFEH.Lbool*VAR.DISP(:) ;
    
end


% Therefore
%
%METHOD = 1;

%
% if METHOD == 1
%
%
%
%
%
%
% end

