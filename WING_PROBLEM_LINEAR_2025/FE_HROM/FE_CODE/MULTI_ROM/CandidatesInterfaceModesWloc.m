function [BasisRdefROT,Vrb_Morth,BasisREFERENCE,TEXTRM,BasisUdomROT] = ...
    CandidatesInterfaceModesWloc(TEXTRM,iface1,iface2,BasisRdef,f1,f2,SingVal_Udef_L,BasisUdef,DATAIN,Vrb,M,BasisUrb)


TEXTRM{end+1} = ['--------------------------------------------'] ;
TEXTRM{end} = ['CANDIDATES FOR BEING INTERFACE MODES']  ;
TEXTRM{end+1} = ['faces = ',num2str(iface1),' and ',num2str(iface2)]  ;


BasisUdef_L1 = BasisUdef(f1,:) ;
BasisUdef_L2 = BasisUdef(f2,:) ;

SingVal_Udef_L = SingVal_Udef_L/SingVal_Udef_L(1) ;

BasisUdef_L1 = bsxfun(@times,BasisUdef_L1',SingVal_Udef_L)' ;
BasisUdef_L2 = bsxfun(@times,BasisUdef_L2',SingVal_Udef_L)' ;

DATAIN = DefaultField(DATAIN,'ROTATION_LOC',[]) ;

DATAIN = DefaultField(DATAIN,'RotationMatrixLocal',[]); 
 
if ~isempty(DATAIN.RotationMatrixLocal)
ROT{1} = DATAIN.RotationMatrixLocal{iface1} ;  % Maps vectors expressed in interface. system to domain system
ROT{2} = DATAIN.RotationMatrixLocal{iface2} ;
else
    ROT{1} = []; ROT{2} = [] ; 
end
 


BasisRdefROT = {BasisRdef(f1,:); BasisRdef(f2,:)} ;
BasisUdomROT = {[BasisUrb(f1,:),BasisUdef(f1,:)], [BasisUrb(f2,:),BasisUdef(f2,:)]} ;
BasisREFERENCE = {BasisUdef_L1, BasisUdef_L2};

% Convert BasisRdefROT,BasisUdomROT and BasisREFERENCE to interface
% coordinates
for iii = 1:length(ROT)
    if ~isempty(ROT{iii}) 
    BasisRdefROT{iii} = RotateMatrix(ROT{iii}',BasisRdefROT{iii})   ; 
    BasisUdomROT{iii} = RotateMatrix(ROT{iii}',BasisUdomROT{iii})   ; 
    BasisREFERENCE{iii} = RotateMatrix(ROT{iii}',BasisREFERENCE{iii})   ; 
    end
end
 BasisREFERENCE = cell2mat(BasisREFERENCE) ; 
    
     

BeamModesInterface = Vrb;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we make the candidates for interface modes M-orthogonal to
% matrix BeamModesInterface
Mchol = chol(M) ;
DATAIN = DefaultField(DATAIN,'TOLERANCE_TRUNCATE_CANDIDATE_INTERFACE_MODES',1e-6) ;  % Most reliable parameter
% for truncating the candidates. 9-Jan-2019
BasisREFERENCE = MorthogonalMatrix(M,Mchol,BasisREFERENCE,BeamModesInterface,DATAIN)  ;
%end
%%%%%%%%%%%%%%%


TEXTRM{end+1} = ['********************************************'] ;
TEXTRM{end+1} = ['Total number of candidates for INTEFACE MODES =  ',num2str(size(Vrb,2)),' + ',num2str(size(BasisREFERENCE,2)),' = ',...
    num2str(size(Vrb,2) + size(BasisREFERENCE,2))] ;

% Now we turn BeamModesInterface (rigid body) M-orthogonal
Xbar = Mchol*BeamModesInterface ;
TOL = 0;
[Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
Vrb_Morth = Mchol\Ubar ;
TEXTRM{end+1} = ['--------------------------------------------'] ;

