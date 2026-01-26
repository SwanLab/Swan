function    [Reactions,WEIGHT ]=   ReactionForces_fromHROM_v2(SNAPvel_loc,SNAPacel_loc,SNAPpk1stress_loc,DISP_LOC,MESH,OPERFE,BasisU,DOFlHROM,...
    DOFrHROM,BasisU_r_proj,BasisPone,ECMdata,DISP_CONDITIONS,Fbody,Ftrac)

if nargin == 0
    load('tmp.mat')
end

% Reconstruction of nodal reactions from HROM variables

% %%%
%  \begin{equation}
%  \label{eq:sdfas--**}
%    \Res{}{\DOFr}  = (\M_{\DOFr} \BasisUall{}{}) \qALLdd{n+1} + (\Ddamp_{\DOFr}\BasisUall{}{}) \qALLd{n+1}  +   \Fint{}{\DOFr}(\qALL{n+1})   -    \Fext{}{\DOFr}(t_{n+1})
%  \end{equation}

%1 ) \M_{\DOFr} \BasisUall{}{}

% 1) INERTIAL FORCES
% ------------------
% ) M_r*BasisUall =  [M_rl  M_rr] diag(BasisU,BasisU_r)
DOFr = DISP_CONDITIONS.DOFr; DOFl = DISP_CONDITIONS.DOFl ;

if ~isempty(SNAPvel_loc)
    M_r = [OPERFE.M(DOFr,DOFl)*BasisU,OPERFE.M(DOFr,DOFr)*BasisU_r_proj   ] ;
    Finer = M_r*SNAPacel_loc ;
else
    Finer = 0 ;
    
end

% 2) DAMPING FORCES
% - -------------------------------
if ~isempty(OPERFE.Ddamp) && ~isempty(SNAPvel_loc)
    D_r = [OPERFE.Ddamp(DOFr,DOFl)*BasisU,OPERFE.Ddamp(DOFr,DOFr)*BasisU_r_proj   ] ;
    Fdamp = D_r*SNAPvel_loc ;
else
    Fdamp =0 ;
end

% 3) INTERNAL FORCES
% ---------------------
% \begin{equation}
%  \Fint{}{\DOFr} =   \overbrace{\Bst{T}{\DOFr} \wST{} \BasisF{}{} \BasisF{+}{\setPoints} }^{\ReconsFintDOFr } \PoneST{}_{\setPoints}   = \ReconsFintDOFr  \PoneST{}_{\setPoints}

% \Bst{T}{\DOFr} \wST{} \BasisF{}{} \BasisF{+}{\setPoints}

% \end{equation}
DATA = [] ;
DATA= DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0);

DATA.MESH.ndim = size(MESH.COOR,2) ;
OPERFE.Bst = OPERFE.Bst(:,DOFr) ;
FintR = cell(1,size(BasisPone,2)) ;
for imodes = 1:size(BasisPone,2)
    FintR{imodes} = InternalForces(OPERFE,BasisPone(:,imodes),[],DATA) ;
end
FintR = cell2mat(FintR)   ;

% RECONSTRUCTION PK1 STRESSES
% -----------------------------

setIndices = small2large(ECMdata.setPoints,DATA.MESH.ndim^2) ;
BasisPoneZ = BasisPone(setIndices,:) ;
coeff =  (BasisPoneZ'*BasisPoneZ)\(BasisPoneZ'*SNAPpk1stress_loc) ;
FintR =  FintR*coeff ;


% External forces
% ---------------
% Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;
FextR =  Fbody.U*Fbody.a +  Ftrac.U*Ftrac.a;

Reactions = Finer + Fdamp + FintR -FextR ;







