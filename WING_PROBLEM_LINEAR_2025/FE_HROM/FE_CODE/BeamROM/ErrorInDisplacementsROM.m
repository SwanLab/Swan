function DIMENSIONLESS_ERROR = ...
    ErrorInDisplacementsROM(DATAROM,MESH1D,qDEF,qRB,DATA_REFMESH,DATAIN,ndimINTF,DOFr,a) ;

% See implementation in
%/home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/
%COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/PAPER_multiLEARN/MeasuringError.tex
% --------------------------------------------------------------------------------------------------

if nargin ==0
    load('tmp.mat')
end

% Definition of matrices U{1} and U{2} (according to the definitions given in MeasuringError.tex)
ientity = 1; 
DATAROM = DATAROM{ientity}; 
DATA_REFMESH = DATA_REFMESH{ientity}; 
% First of all, the basis matrix of displacements 
BasisUdom = [DATA_REFMESH.BasisUrb,DATAROM.BasisUdef]; 
% Now the rotation matrices 
% ........................
ROT = DATA_REFMESH.RotationMatrixFace  ; 
f = {DATAROM.f1,DATAROM.f2} ; 
%BasisUdomF = cell(size(f)) ; 
U = cell(size(f)) ; 
W = U ; % For boundaries 

ndim  = size(DATA_REFMESH.COOR,2) ; 
% Mass matrix interfaces 
M1d = DATA_REFMESH.GeometricMassMatrixInterface ;
% We set the same for both 
M1d{2} = M1d{1} ; 
Mintf = cell(length(f),1) ;
Mchol = cell(length(f),1) ; 
for iface= 1:length(f)
    Mintf{iface} = sparse(size(M1d{iface},1)*ndim,size(M1d{iface},2)*ndim) ;
    for idim = 1:ndim
        Mintf{iface}(idim:ndim:end,idim:ndim:end) = M1d{iface} ;
    end
    Mchol{iface} = chol( Mintf{iface}) ; 
    
end


for iface = 1:length(f)
    if isempty(ROT{iface})
        BasisUdomF  = BasisUdom(f{iface},:) ; 
    else
        % Rotation matrix ------------------------------------------------
        BasisUdomF  = RotateMatrix(ROT{iface},BasisUdom(f{iface},:)) ;        
    end
    % Mass matrix of the interface
     U{iface} = Mchol{iface}*BasisUdomF ; 
     W{iface} = Mchol{iface}*DATAROM.BasisInt ; 
end

% Now we compute   
UQUAD=  {U{1}'*U{1} , -U{1}'*U{2}
         -U{2}'*U{1} , U{2}'*U{2}}  ; 
     
UWQUAD11 =       {U{1}'*U{1} , -U{1}'*W{2}
         -W{2}'*U{1} , W{2}'*W{2}}  ; 

UWQUAD22 =       {W{1}'*W{1} , -W{1}'*U{2}
         -U{2}'*W{1} , U{2}'*U{2}}  ;      
     
     
 AREA = UWQUAD22{1,1}(1) ;  
     
 % Loop over  1D nodes 
 
 ERROR2 = zeros(size(MESH1D.COOR,1),1) ;
 DISP2 = zeros(size(MESH1D.COOR,1),1) ;
 
 for  inode = 1:size(MESH1D.COOR,1)
     [ELEMS,INODEloc] = find(MESH1D.CN==inode) ; 
     if length(ELEMS)==1 
         % Boundary node 
         % Associated DOFs 
         DOFloc= small2large(inode,ndimINTF) ; 
         % Check if it is prescribed
         III = intersect(DOFr,DOFloc) ; 
         if length(III) == length(DOFloc)
             aLOC = a(DOFloc) ;
             qDOM = cell(2,1) ;
             qDOM{INODEloc} = [qRB(:,ELEMS); qDEF(:,ELEMS)] ; 
             IOTHER = setdiff(1:2,INODEloc) ; 
             qDOM{IOTHER} = aLOC ; 
             if INODEloc == 1
                 UU = UWQUAD11 ; 
                  DISP2(inode) =  qDOM{2}'*UU{2,2}*qDOM{2}   ; 
             else
                 UU = UWQUAD22 ; 
                 DISP2(inode) =  qDOM{1}'*UU{1,1}*qDOM{1}   ; 
             end
             %%% 
             ERROR2(inode) = cell2mat(qDOM)'*cell2mat(UU)*cell2mat(qDOM) ; 
           
             
         else
             % Error cannot be computed. Free node (or partially free)
             ERROR2(inode) = 0 ; 
            DISP2(inode) = 0; % REference 
         end
     else
         % INTERIOR NODE 
         % -------------
         qDOM = cell(2,1) ; 
         for idom  =1:2
             qDOM{INODEloc(idom)} =  [qRB(:,ELEMS(idom)) ; qDEF(:,ELEMS(idom))  ] ; 
         end
         UU = UQUAD ; 
         
         ERROR2(inode) = cell2mat(qDOM)'*cell2mat(UU)*cell2mat(qDOM) ; 
         DISP2(inode) = 0.5*( qDOM{1}'*UU{1,1}*qDOM{1} + qDOM{2}'*UU{2,2}*qDOM{2})  ; 
         
         
         
     end
     
       
     
     
 end 
 
  
%  % DIMENSIONLESS -> Divided by the maximun displacement 
%  a = reshape(a,ndimINTF,[]) ; 
%  aTRANS = sqrt(sum(a(1:3,:).^2,1)) ; 
%  aMAX = max(aTRANS) ; 

DATAIN = DefaultField(DATAIN,'ShowErrorDisplacementsInPercentage',1) ; 
 
if DATAIN.ShowErrorDisplacementsInPercentage == 1
DIMENSIONLESS_ERROR = sqrt(ERROR2/max(DISP2))*100 ;
else
    DIMENSIONLESS_ERROR = sqrt(ERROR2/AREA)*1000 ; 
end

 %



 



end
