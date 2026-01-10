function [V,DATAOUT] = OutputKinemApproach(Dcomp,DOFs,DOFm,BasisINTFall,nDOFsFACEall,DATAOUT,IndicesRB,BasisINTFcand)


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