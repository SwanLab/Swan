function  [DATAOUT,V,TEXTP] =  InterfaceRVEenergetic(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
    SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,Mdom,SingVal_Udef)


% STEP 1
% -------
% Determination of BasisINTFall along with the indexes of
% master/slaves,rigid-body DOFS ---- (candidates)
% Matrix Dcomp--- 
MasterDOFS_perface = {1e10,1e10} ; % First iteration, this is set to a very high value 

 [BasisINTFall,Indexes,TEXTP,BasisINTFcand,IndicesRB,Dcomp,nDOFsFACEall] ...
     = BasisINTFall_local(Vrb,FACES_GROUPS,BasisUdef,BasisRdef,M,DATAIN,...
    SingVal_Udef,SinvVal_Rdef,TEXTP,MasterDOFS_perface,fI,BasisUrb,Mdom) ; 

% ------------------------------------------------------------------------ 
% STEP 2. Master and Slave DOFs
% Select linearly independent rows. 
% --------------------------------------

 [DOFm,DOFs] = CoarseDOFS_DEIMbased(Dcomp,DATAIN,IndicesRB,nDOFsFACEall) ;

%[DOFm,DOFs] = SelectDOFsMASTER_SLAVES_methodENER(Dcomp,Indexes) ; 


%[DOFm,DOFs] = CoarseDOFS_master_slave(Dcomp,DATAIN,BasisUdom_f,IndicesRB,IndicesDEF,nDOFsFACEall,nRB) ;

TEXTP{end+1} = ['MASTER DOFs = ',num2str(DOFm(:)')] ;

DOFmRB = intersect(IndicesRB,DOFm) ;
TEXTP{end+1} = ['MASTER DOFs (of rigid body type) = ',num2str(DOFmRB(:)')] ;

if rank(Dcomp(DOFm,:)) ~=length(DOFm)
    error('The selection of Master DOFs is conducive to an ill-conditioned Dcomp(DOFm,:)')
end

Acomp = Dcomp(DOFs,:)*inv(Dcomp(DOFm,:)) ;

% OUTPUT WILL BE A STRUCTURED ARRAY

V.DOFmP = DOFm ;   % Master DOFs
V.DOFsP = DOFs ;   % Slave DOFs
V.IndicesRB = IndicesRB ; % Indexes Rigid Body DOFs
V.BasisINTFall_cell = BasisINTFcand ;  % Basis containing all modes (cell array)
V.BasisINTF = BasisINTFall(:,DOFm) + BasisINTFall(:,DOFs)*Acomp ;   % Basis matrix to be used in the formulatiion
V.Acomp = Acomp ;
% Number of  DOFs per face
iacum = 1;
DOFsFACE = cell(size(nDOFsFACEall)) ;
DOFsFACE_V = cell(size(nDOFsFACEall)) ;
iacum_V = 1;
for iface = 1:length(nDOFsFACEall)
    DOFsFACEloc = iacum:(iacum+nDOFsFACEall(iface)-1) ;
    DOFsFACE{iface} = intersect(DOFsFACEloc,DOFm) ;
    DOFsFACE_V{iface} =  iacum_V:(iacum_V+length(DOFsFACE{iface})-1) ; ;
    iacum = DOFsFACEloc(end)+1 ;
    iacum_V = DOFsFACE_V{iface}(end)+1 ;
end
V.DOFsFACE = DOFsFACE ;
V.DOFsFACE_V = DOFsFACE_V ;

[nDOFsFACE ]= cellfun(@length,DOFsFACE);
V.nDOFsFACE = nDOFsFACE ;
DATAOUT.BasisInt = V  ;


