function [FextPR, tPRESSst,kPRESSst ]= PressureForcesFromDisplacements(DATA,OPERFE,VAR) 
% This function returns the Gauss points forces due to hydrostatic pressure
% (tPRESSst), the contribution to the tangent matrix to such forces
% (kPRESSst), 
% and the corresponding nodal force vector FextPR 
% See formulation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
%  JAHO, 1-July-2021, Cartagena, SPAIN (UPCT)

if nargin == 0
    load('tmp.mat')
                    DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO = 1 ; 

end



% %  \begin{equation}
%   \label{eq:83osdd}
%    \tPRESS^e_{\xiFEft_g}  =     (p \mNORM)_{\xiFEft_g}
%   \end{equation}

% 1. Normal at each gauss point (not unitary), and contribution to the
% stiffness matrix 
% -------------------------------------------
[mNORMst,kNORMst] = NormalNonUnitaryBoundary(OPERFE.HYDRO,VAR,DATA)  ; 

% 2. Hydrostatic pressure 
%--------------------
[pST,kPst,Y_neg]=  HydrostaticPressure(OPERFE.HYDRO,VAR,DATA) ; 

if  DATA.FOLLOWER_LOADS.HYDROSTATIC.INCLUDE_DYNAMIC_PRESSURE == 1
    [pDYN,kPdyn,vST]=  DynamicPressure(OPERFE.HYDRO,VAR,DATA,mNORMst,kNORMst,Y_neg) ;
    pST = pST + pDYN ;
    kPst = kPst + kPdyn ;
else
    ndim = 3;
    vST = ConvertBlockDiag_general(OPERFE.HYDRO.NbST,ndim,OPERFE.HYDRO.irowsNDIM,OPERFE.HYDRO.icolsNDIM)*OPERFE.HYDRO.Lbool*VAR.VEL(:) ;
    
end




ndim = length(mNORMst)/length(pST) ; 

% Force 
% -------
tPRESSst = zeros(size(mNORMst)) ; 

for idim = 1:ndim 
    tPRESSst(idim:ndim:end) = pST.*mNORMst(idim:ndim:end)  ; 
end 


if  DATA.FOLLOWER_LOADS.HYDROSTATIC.WET_DAMPING_COEFFICIENT ~= 0    
    % Wet damping  == Standard damping , but only applied to wet surfaces 
    % -------------------------------------------------------------------
   % t_wet = mu_wet*v*HEAVISIDE(y)
    mu_wet = DATA.FOLLOWER_LOADS.HYDROSTATIC.WET_DAMPING_COEFFICIENT ; 
    
    tWETdamping = zeros(size(mNORMst)) ;  
    BIN_yNEG= Y_neg < 0 ; 
    
    for idim = 1:ndim 
        tWETdamping(idim:ndim:end) = mu_wet*vST(idim:ndim:end).*BIN_yNEG.*OPERFE.HYDRO.JacobianWeights(:); 
     %   tWETdamping(idim:ndim:end) = mu_wet*vST(idim:ndim:end).*OPERFE.HYDRO.JacobianWeights(:); 
    end
        FextPR = -OPERFE.HYDRO.NbST_w'*(tPRESSst+tWETdamping); 

else
    
    FextPR = -OPERFE.HYDRO.NbST_w'*(tPRESSst); 

    
end




% RESULTING NODAL FORCES 
% 
% \begin{equation}
% \label{eq:d4555}
%  \boxed{\FextPR{}{}   = - \NbST{T}  \wSTft{}  \tPRESSst}
% \end{equation}



% Contribution to tangent matrix 
% ------------------------------
% \begin{equation}
%   [\T^e]_i  = \overbrace{[\mNORMst^e]_i   \hadmP  \derpar{\pPRESSst^e}{\d^e}}^{\ddSS [\W^e]_{i}}   +  \overbrace{\pPRESSst^e \hadmP [\derpar{\mNORMst^e}{\d^e}]}^{\ddSS [\V^e]_i}
% \end{equation}

if   DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO  == 0
    
    kPRESSst = zeros(size(kNORMst))  ;
    
    for idim = 1:ndim
        kPRESSst(idim:ndim:end,:) =   bsxfun(@times,kNORMst(idim:ndim:end,:), pST) ;
        kPRESSst(idim:ndim:end,:)  = kPRESSst(idim:ndim:end,:)  + bsxfun(@times,kPst, mNORMst(idim:ndim:end)) ;
        
        if   DATA.FOLLOWER_LOADS.HYDROSTATIC.WET_DAMPING_COEFFICIENT ~= 0
            
            kPRESSst(idim:ndim:end,:) =  kPRESSst(idim:ndim:end,:) + mu_wet*bsxfun(@times,OPERFE.HYDRO.NbST(idim:ndim:end,:),BIN_yNEG.*OPERFE.HYDRO.JacobianWeights(:)) ;
            
            %     kPRESSst(idim:ndim:end,:) =  kPRESSst(idim:ndim:end,:) + mu_wet*bsxfun(@times,OPERFE.HYDRO.NbST(idim:ndim:end,:),OPERFE.HYDRO.JacobianWeights(:)) ;
            
        end
    end
    
else
    kPRESSst = [ ] ; 
end



 


