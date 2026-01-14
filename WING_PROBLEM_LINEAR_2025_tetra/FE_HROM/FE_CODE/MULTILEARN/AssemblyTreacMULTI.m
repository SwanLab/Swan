
function  [T, c, s]=  AssemblyTreacMULTI(BasisRdef,f1,f2,nDOM,reactDOMrb,ALPHA,COL,COLloc)

%dbstop('4')
if nargin == 0
    load('tmp3.mat')
end


f = [f1 ; f2] ;
% Cov_f = BasisRdef(f,:)'*BasisRdef(f,:) ;
% Cov_f2f1 = BasisRdef(f2,:)'*BasisRdef(f1,:) ;
% Cov_f1f2 = BasisRdef(f1,:)'*BasisRdef(f2,:) ;
% Cov_f2f2 = BasisRdef(f2,:)'*BasisRdef(f2,:) ;
% Cov_f1f1 = BasisRdef(f1,:)'*BasisRdef(f1,:) ;

%nMODES = size(BasisRdef,2) ;
%nrows = nMODES*nDOM ;
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nMODES ]= cellfun(@size,BasisRdef) ;  % Number of deformation modes
%nMODESrve = nDEFrve + nRB ; % Total number of modes
nrows  = sum(nDOMrve.*nMODES) ;
T = sparse(nrows,nrows) ;
c = zeros(nrows,1);
s = 0 ;
iROWS = 0 ;
iCOLS = 0 ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) ;
end
for idom = 1:nDOM
    itype_i = COLloc(idom) ;  % Type of RVE
    Cov_f1f1 = BasisRdef{itype_i}(f1,:)'*BasisRdef{itype_i}(f1,:) ;
    Cov_f2f2 = BasisRdef{itype_i}(f2,:)'*BasisRdef{itype_i}(f2,:) ;
    Cov_f = BasisRdef{itype_i}(f,:)'*BasisRdef{itype_i}(f,:) ;
    if idom == 1
        if nDOM >1
            
            itype_ip1 = COLloc(idom+1) ;
            Cov_f2f1 = BasisRdef{itype_i}(f2,:)'*BasisRdef{itype_ip1}(f1,:) ;
            VECT = [(1-ALPHA(1))*Cov_f1f1 +  Cov_f2f2, Cov_f2f1] ;
            indROWS = 1:size(VECT,1) ;
            indCOLS = 1:size(VECT,2) ;
            T(indROWS,indCOLS) = VECT ;
            c(indROWS) = (1-ALPHA(1))*(reactDOMrb{idom}(f1,:))'*BasisRdef{itype_i}(f1,:)  + ...
                (reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))'*BasisRdef{itype_i}(f2,:) ;
            s = s + norm(reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))^2  ...
                +  (1-ALPHA(1))*norm(reactDOMrb{idom}(f1,:))^2;
            iROWS = length(indROWS) ; iCOLS = 0 ; 
        else
            VECT=  [(1-ALPHA(1))*Cov_f1f1 + (1-ALPHA(end))*Cov_f2f2] ; 
            indROWS = 1:size(VECT,1) ;
            indCOLS = 1:size(VECT,2) ;
            T(indROWS,indCOLS) = VECT ; 
            c(indROWS) =(1-ALPHA(1))*(reactDOMrb{idom}(f1,:) )'*BasisRdef(f1,:) + ...
                (1-ALPHA(end))*(reactDOMrb{idom}(f2,:) )'*BasisRdef(f2,:) ;
            s = s + (1-ALPHA(1))*norm(reactDOMrb{idom}(f1,:))^2 + ...
                (1-ALPHA(end))*norm(reactDOMrb{idom}(f2,:))^2  ;
        end
    elseif idom == nDOM
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        Cov_f1f2 = BasisRdef{itype_i}(f1,:)'*BasisRdef{itype_im1}(f2,:) ;%
        VECT = [Cov_f1f2, Cov_f1f1+ (1-ALPHA(end))*Cov_f2f2] ; ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        T(indROWS,indCOLS) =  VECT ; 
        c(indROWS) = (reactDOMrb{idom-1}(f2,:) + reactDOMrb{idom}(f1,:))'*BasisRdef{itype_i}(f1,:) + ...
            (1-ALPHA(end))*reactDOMrb{idom}(f2,:)'*BasisRdef{itype_i}(f2,:);
        s = s + (1-ALPHA(end))*norm(reactDOMrb{idom}(f2,:))^2 ;
    else
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        itype_ip1 = COLloc(idom+1) ; % Type of contiguous RVE
        Cov_f2f1 = BasisRdef{itype_i}(f2,:)'*BasisRdef{itype_ip1}(f1,:) ;
        Cov_f1f2 = BasisRdef{itype_i}(f1,:)'*BasisRdef{itype_im1}(f2,:) ;%
        VECT = [Cov_f1f2, Cov_f, Cov_f2f1] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        T(indROWS,indCOLS) = [Cov_f1f2, Cov_f, Cov_f2f1] ;
        c(indROWS) = (reactDOMrb{idom-1}(f2,:) + reactDOMrb{idom}(f1,:))'*BasisRdef{itype_i}(f1,:) + ...
            (reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))'*BasisRdef{itype_i}(f2,:)   ;
        s = s + norm(reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))^2 ;
        iROWS = iROWS + length(indROWS) ;
        iCOLS = sum(nMODES_all(1:idom-1)) ;
        
    end
end



end





