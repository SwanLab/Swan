function [BasisINT,TEXTR] =  ChooseIntModesWORK1(Vrb_Morth,Vdef_Morth,SingVal_Rdef,BasisRdefROT,TEXTR,nmodesR,Vrb)


U  = [Vrb_Morth,Vdef_Morth ];  % Candidates for being interface modes
% Reaction force matrix (weighted by the SINGULAR VALUES )
SingVal_Rdef = SingVal_Rdef/SingVal_Rdef(1) ;
% Rotated reactions weighted by the singular values
nfaces =2 ;
rancTall = [] ; 
for i = 1:nfaces
    % Matrix T for all candidates (weighted by the singular values)    
    BasisRdefROT_W{i} = bsxfun(@times,BasisRdefROT{i}',SingVal_Rdef)' ;
    CW{i} = BasisRdefROT_W{i}'*U ;    
    C = BasisRdefROT{i}'*U ;    
    TOL  = 1e-3 ; 
    DATLOC.RELATIVE_SVD = 1; 
    [~,SSS,~] = SVDT(C,TOL,DATLOC ) ; 
    rancTall(i) = length(SSS) ;  %
end

nmodesR = min(nmodesR,min(rancTall)) ; 

% The problem boils down to select independent rows of C{i} i=1,2
Irb = 1:size(Vrb_Morth,2) ; % Indices RIGID BODY
% Check rank RIGID BODY MODES
% ----------------------------
for i = 1:nfaces
    CrbLOC = BasisRdefROT{i}'*U(:,Irb) ;
    if rank(CrbLOC) ~= length(Irb)
        error('Not proper training ---rigid body modes cannot be transmittedb')
    end    
end
%%%%
% Initialization
TOL_ERROR_APPROX  = 1e20 ;
k = length(Irb) ;
I = Irb ;
while    k < nmodesR
    nrW = cell(size(CW)) ;
    for i = 1:nfaces
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
TEXTR{end+1} = '***************************' ; 
I = sort(I) ;
TEXTR{end+1} = ['SELECTED INTERFACE MODES   = ',num2str(I)] ;
TEXTR{end+1} = '***************************' ; 
Idef = I((length(Irb)+1):end) ;
BasisINT = [Vrb,U(:,Idef)] ;
