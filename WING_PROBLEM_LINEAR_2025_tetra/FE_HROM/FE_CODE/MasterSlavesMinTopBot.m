function [MASTER_uvFL,SLAVE_uvFL,G_uvFL,uBmacro] = ...
    MasterSlavesMinTopBot(NODESpl,CNb,COOR,TypeElementB,strainINP,DATA,h,NODESzBOT,NODESzTOP)
% See Plates_comp_homogenization.pdf
%dbstop('5')
if nargin==0
    load('tmp.mat')
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
% Reference node (SLAVE)
iREF =1; 
nodzREF = nodesZtop(iREF) ; 
if abs(Atop(iREF))<1e-16
    error('Choose another reference/slave point')
end
AtopREF = Atop(iREF) ; 
% Thus 
nodesREM = nodesZtop ; 
nodesREM(iREF) = [] ; 
MASTER_uvFL = [nodesREM; nodesZbot] ; 
SLAVE_uvFL = nodzREF ; 

% Matrix "G" corresponding to these master/slave nodes
REMnodes =1:length(nodesZtop) ; 
REMnodes(iREF) = []; 
nmastUV = length(MASTER_uvFL) ; 
nslaveUV = 1; 
G_uvFL = [-Atop(REMnodes) Abot]/AtopREF ; 

% COORslv = COOR(SLAVE_uvFL,:) ; 
% COORslv = reshape(COORslv,3,[]) ; 
% COORmst = COOR(MASTER_uvFL,:) ; 
% COORmst = reshape(COORmst,3,[]) ; 
% 
% COORmin.MASTER = COORmst ; 
% COORmin.SLAVE =  COORslv;

%%%%
COORtop = COOR(nodesZtop,:) ; 
%COORtop = reshape(COORtop',3,[]) ; 
dB_top= CoarseDisplace(strainINP,COORtop',DATA,h,[]) ;
COORbot = COOR(nodesZbot,:) ; 
%COORbot = reshape(COORbot',3,[]) ; 
dB_bot= CoarseDisplace(strainINP,COORbot',DATA,h,[]) ;

%% Therefore, for u
idim=1 ; 
uBmacro.u = (Atop*dB_top(idim,:)' - Abot*dB_bot(idim,:)')/AtopREF ;
idim=2 ; 
uBmacro.v = (Atop*dB_top(idim,:)' - Abot*dB_bot(idim,:)')/AtopREF ;


 