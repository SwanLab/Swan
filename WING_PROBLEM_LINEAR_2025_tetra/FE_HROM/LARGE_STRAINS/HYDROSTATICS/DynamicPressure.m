function [pDYN,kPdyn,vST]=  DynamicPressure(OPERHYDRO,VAR,DATA,mNORMst,kNORMst,Y_neg) 
% % % Dynamic pressure (waterline X2 = 0).  pDYN =
% 0.5*densW*v_n^2*HEAV(v_n)*Y_neg 
% % % See formulation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
% %  JAHO, 3-July-2021, Cartagena, SPAIN (UPCT)
% ----------------------------------------------------------------------------------------------------------------------
if nargin == 0 
    load('tmp.mat')
end

%  \begin{equation}
%  \v_{\xiV} = (\n^T  \NlocBm{e}{}_{\xiV})   \dD^e  
% \end{equation}

 densW = DATA.FOLLOWER_LOADS.HYDROSTATIC.densW ; 

% 1. UNIT NORMAL  
% --------------
ndim  = DATA.MESH.ndim ; 
nST = reshape(mNORMst,ndim,[]) ;  
nST_norm = sqrt(sum(nST.^2,1)) ;  % Norm 
nST = bsxfun(@times,nST',1./nST_norm')'; 
nST = nST(:); 
% ---------------------------------------------------


% 2. Velocities at all the Gauss points of the boundary
% ----------------------------------
% vST = Nshape*Lbool*VELOCITY_NODES
vST = ConvertBlockDiag_general(OPERHYDRO.NbST,ndim,OPERHYDRO.irowsNDIM,OPERHYDRO.icolsNDIM)*OPERHYDRO.Lbool*VAR.VEL(:) ; 


% 3. Normal component of the velocity 
% ------------------------------------
vNst = zeros(size(nST_norm')) ;  
for idim = 1:ndim 
    vNst =  vNst + nST(idim:ndim:end).*vST(idim:ndim:end)  ; 
end

% 4. Set to zero all entries whose component along the normal are negative 
% -------------------------------------------------
vNst_pos = vNst.*(vNst>0) ; 


% 5. Dynamic pressure 
% -------------------
% 0.5*densW*v_n^2*HEAV(v_n)*Y_neg 

    BIN_yNEG= Y_neg < 0 ; 


% vNst_pos = vNst_pos.*Y_neg ; 

 vNst_pos = vNst_pos.*BIN_yNEG ; 

pDYN = 0.5*densW*(vNst_pos.^2)  ; 


% ***************************************************************+
%*******************************************************************

% CONTRIBUTION TO STIFFNESS TANGENT MATRIX 
% ------------------------------------------

%   \begin{equation}
%  \label{eq:4.m}
%  \k_{dyn} =    \densW  v_n  \cblue{\derpar{v_n}{\d^e} }\HEAV{v_n} \HEAV{y}
%  \end{equation} 

% Let's start by \derpar{v_n}{\d^e}  
if  DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO == 0
kPdyn = StiffnessDynamicPressure(kNORMst,nST_norm',mNORMst,DATA,vST,vNst_pos,nST,OPERHYDRO,densW) ; 
else
    kPdyn = 0 ; 
end




 
