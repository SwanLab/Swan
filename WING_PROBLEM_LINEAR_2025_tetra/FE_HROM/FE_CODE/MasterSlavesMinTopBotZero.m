function [MASTER_uvFLglo,SLAVE_uvFLglo,G_uvFLglo,uBmacroGLO] = ...
    MasterSlavesMinTopBotZero(NODESpl,CNb,COOR,TypeElementB,strainINP,DATA,h,NODESzBOT,NODESzTOP)
% See Plates_comp_homogenization.pdf
%dbstop('5')
if nargin==0
    load('tmp1.mat')
end
%% Matrices Atop and Abot
% --------------------------
% Atop 
% -----
% Nodes pertaining to plane  z = zmax , 

nodesZtop = NODESpl{3,1} ; % Excluding lines 
 %nodesZtop = NODESzTOP ; % Including lines 
Atop = AtopBot_cal(CNb,nodesZtop,COOR,TypeElementB,NODESzTOP) ; 
% Abot
% -----
% Nodes pertaining to plane  z = zmin ,
nodesZbot = NODESpl{3,2} ; 
 %nodesZtop = NODESzBOT ; % Including lines 
Abot = AtopBot_cal(CNb,nodesZbot,COOR,TypeElementB,NODESzBOT) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference node (SLAVE) TOP 
% ---------------------------
iREF =1; 
nodzREF = nodesZtop(iREF) ; 
if abs(Atop(iREF))<1e-16
    error('Choose another reference/slave point')
end
AtopREF = Atop(iREF) ; 
% Thus 
nodesREM = nodesZtop ; 
nodesREM(iREF) = [] ; 
MASTER_uvFL.top = [nodesREM] ; 
SLAVE_uvFL.top = nodzREF ; 
% Reference node (SLAVE) BOTTOM 
% ---------------------------
iREF =1; 
nodzREF = nodesZbot(iREF) ; 
if abs(Abot(iREF))<1e-16
    error('Choose another reference/slave point')
end
AbotREF = Abot(iREF) ; 
% Thus 
nodesREM = nodesZbot ; 
nodesREM(iREF) = [] ; 
MASTER_uvFL.bot = [nodesREM] ; 
SLAVE_uvFL.bot = nodzREF ; 


% Matrix "G" corresponding to these master/slave nodes 
% TOP 
% --------------------------------------------------------
REMnodes =1:length(nodesZtop) ; 
REMnodes(iREF) = []; 
nmastUV.top = length(MASTER_uvFL.top) ; 
nslaveUV.top = 1; 
G_uvFL.top = [-Atop(REMnodes) ]/AtopREF ; 
% BOTTOM
% --------------------------------------------------------
REMnodes =1:length(nodesZbot) ; 
REMnodes(iREF) = []; 
nmastUV.bot = length(MASTER_uvFL.bot) ; 
nslaveUV.bot = 1; 
G_uvFL.bot = [-Abot(REMnodes) ]/AbotREF ; 

%%%%
COORtop = COOR(nodesZtop,:) ; 
COORtop = reshape(COORtop',3,[]) ; 
dB_top= CoarseDisplace(strainINP,COORtop,DATA,h,[]) ;
COORbot = COOR(nodesZbot,:) ; 
COORbot = reshape(COORbot',3,[]) ; 
dB_bot= CoarseDisplace(strainINP,COORbot,DATA,h,[]) ;

%% Therefore, for u
% TOP 
idim=1 ; 
uBmacro.u.top = (Atop*dB_top(idim,:)')/AtopREF ;
idim=2 ; 
uBmacro.v.top = (Atop*dB_top(idim,:)' )/AtopREF ;
% BOTTOM
idim=1 ; 
uBmacro.u.bot = (Abot*dB_bot(idim,:)')/AbotREF ;
idim=2 ; 
uBmacro.v.bot = (Abot*dB_bot(idim,:)' )/AbotREF ;


%%%% 
MASTER_uvFLglo = [MASTER_uvFL.top; MASTER_uvFL.bot]; 
SLAVE_uvFLglo = [SLAVE_uvFL.top ; SLAVE_uvFL.bot] ; 
 G_uvFLglo = sparse(length(SLAVE_uvFLglo),length(MASTER_uvFLglo))  ; 
 
G_uvFLglo(1:nslaveUV.top, 1:nmastUV.top) = G_uvFL.top ; 
G_uvFLglo(nslaveUV.top+1:end, nmastUV.top+1:end) = G_uvFL.bot ; 
%%
uBmacroGLO.u = [uBmacro.u.top ; uBmacro.u.bot]  ; 
uBmacroGLO.v = [uBmacro.v.top ; uBmacro.v.bot]  ; 
