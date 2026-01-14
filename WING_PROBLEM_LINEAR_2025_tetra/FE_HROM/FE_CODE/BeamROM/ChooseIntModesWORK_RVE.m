function [BasisINT,TEXTR] =  ChooseIntModesWORK_RVE(BasisINTFall,SingVal_Rdef,...
    BasisRdefROT_all,TEXTR,Vrb,FACES_GROUPS,DATAIN)

% Adaption of ChooseIntModesWORK.m for the case of RVEs (4 faces). 
% 17-May-2019 
if nargin == 0
    load('tmp1.mat')
end
DATALOC = DATAIN.INTERFACE_MODES_REACTIONS_WORK ;

DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_DEF_MODES',1e-6) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_ANGLE_INTERSECTION_REACTIONS',0.01) ; % TOLERANCE FOR SVD DEFORMATIONAL MODES
DATALOC = DefaultField(DATALOC,'TOLERANCE_SVD_REACTION_MODES_INTERSECTION',1e-6) ; % TOLERANCE FOR SVD reaction MODES


% --------------------------------------
% Candidates for being interface modes
% --------------------------------------
nRB = size(Vrb{1},2) ; 
faceM(1) = FACES_GROUPS{1}(1) ;   % Minus faces 
faceM(2) = FACES_GROUPS{2}(1) ; 
faceP(1) = FACES_GROUPS{1}(2) ;   % Plus faces  
faceP(2) = FACES_GROUPS{2}(2) ; 

U = BasisINTFall(faceM) ;  
% Rigid body indexes 
Irb = [] ; 
iacum = 0 ; 
Indface = cell(size(U)) ;   % All indexes of each face
IndDEF =  cell(size(U)) ;   % Deformational indices of each face 
DOFSface = IndDEF ; 
iacumDOFS = 0 ; 
for iface = 1:length(U)
    jjj = 1:size(U{iface},2) ;  
    iii= 1:nRB ; 
    kkk = nRB+1:size(U{iface},2) ; 
    Irb = [Irb,iii+iacum] ; 
    Indface{iface} = [jjj+iacum] ; 
    IndDEF{iface} = [kkk+iacum] ; 
    iacum = iacum+ size(U{iface},2) ; 
    DOFSface{iface} = [1:size(U{iface},1)] +  iacumDOFS ; 
    iacumDOFS = iacumDOFS + size(U{iface},1) ; 
end
% Convert U into a diagonal block matrix 
U = blkdiag(U{:}) ; 



% Rotated reactions 

BasisRdefROT = cell(1,2) ; 
BasisRdefROT{1} = [BasisRdefROT_all{faceM(1)}; BasisRdefROT_all{faceM(2)}] ;
BasisRdefROT{2} = [BasisRdefROT_all{faceP(1)}; BasisRdefROT_all{faceP(2)}] ;

% COMPUTING THE INTERSECTION 
% ---------------------------------
INTERS = 0 ; 
if INTERS ==1
    error('It makes no sense to talk about intersection ....')
TOL_SVD = DATALOC.TOLERANCE_SVD_REACTION_MODES_INTERSECTION;
TOL_ANGLE = DATALOC.TOLERANCE_ANGLE_INTERSECTION_REACTIONS ;
[RR,ANGLES ]= IntersectionSubspaces(BasisRdefROT{1},BasisRdefROT{2},TOL_SVD,TOL_ANGLE) ;
nmodesR = (size(RR,2)) ;

    TEXTR{end+1} = '****************************************' ;

   TEXTR{end+1} = ['Dimension of the intersection of reaction subspaces (two pairs of faces) (for TOL ANGLE = ',num2str(TOL_ANGLE),') --> ',num2str(nmodesR)] ;
    TEXTR{end+1} = '****************************************' ;
else
   nmodesR = 1e20 ;  
end


% Reaction force matrix (weighted by the SINGULAR VALUES )
SingVal_Rdef = SingVal_Rdef/SingVal_Rdef(1) ;

% Rotated reactions weighted by the singular values
RankTmatrices  =[]; 
INDICES = {'1_2','3_4'}; 
TEXTR{end+1} = '***************************';

for i = 1:length(BasisRdefROT)
    % Matrix T for all candidates (weighted by the singular values)    
    BasisRdefROT_W{i} = bsxfun(@times,BasisRdefROT{i}',SingVal_Rdef)' ;
    CW{i} = BasisRdefROT_W{i}'*U ;    
    C = BasisRdefROT{i}'*U ; 
    TOL = 1e-3 ;   % To determine the effective rank of the matrices 
    DATALOC.RELATIVE_SVD = 1; 
    [~,SSS,~] = SVDT(C,TOL,DATALOC) ; 
    RankTmatrices(i) = length(SSS) ; 
    TEXTR{end+1} = ['Rank T_  ',INDICES{i},'=',num2str(RankTmatrices(i))]; 

end

nmodesR = min(nmodesR,size(U,2)) ; % Maximum number of  candidate modes 
TEXTR{end+1} = ['Total number of candidate modes = ',num2str(nmodesR)]; 
nmodesR = min(nmodesR,min(RankTmatrices)) ; 
TEXTR{end+1} = ['Total number of interface modes (two faces)= ',num2str(nmodesR)]; 



 % Check rank RIGID BODY MODES
% ----------------------------
for i = 1:length(BasisRdefROT)
    CrbLOC = BasisRdefROT{i}'*U(:,Irb) ;
    if rank(CrbLOC) ~= length(Irb)
        error('Not proper training ---rigid body modes cannot be transmited to neighboring domains')
    end    
end
%%%%
% Initialization
TOL_ERROR_APPROX  = 1e20 ;
k = length(Irb) ;
I = Irb ;
while    k < nmodesR
    nrW = cell(size(CW)) ;
    for i = 1:length(CW)
        % Residual
        rW  = CW{i} - CW{i}(:,I)*(CW{i}(:,I)\CW{i})  ;
        % Norm of the residual
        nrW{i} = sqrt(sum(rW.^2,1))' ;
    end
    nrMW = cell2mat(nrW) ;
    CRIT = prod(nrMW,2) ;
    [~,Inew] = max(CRIT) ;
    I = [I,Inew] ;
    k = k+1;        
end



%% Interface modes
TEXTR{end+1} = '***************************';
I = sort(I) ;
TEXTR{end+1} = ['SELECTED INTERFACE MODES  (2 faces ) = ',num2str(I)] ;
TEXTR{end+1} = '***************************';

BasisINT = cell(size(Vrb)) ; 
for iface = 1:length(Vrb)
    IM =   find(iface == faceM) ; % Face minus
    IP =   find(iface == faceP) ; % Face plus     
    IND = IM ; 
    if isempty(IM)
        IND = IP ;  % Indices for choosing which face
    end
    IndDEFloc = IndDEF{IND} ;  % Deformational modes of this face (all)
    IndDEFloc_inter = intersect(I,IndDEFloc)  ;   % intersection with the chosen  modes
    BasisINT{iface} = [Vrb{iface},U(DOFSface{IND},IndDEFloc_inter)] ; 
end

 
