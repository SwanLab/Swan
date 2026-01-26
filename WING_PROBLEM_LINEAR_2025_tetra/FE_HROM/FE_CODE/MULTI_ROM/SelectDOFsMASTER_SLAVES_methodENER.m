function [DOFm,DOFs] = SelectDOFsMASTER_SLAVES_methodENER(Dcomp,Indexes) ; 

if nargin ==0 
    load('tmp1.mat')
end

% Total number of candidates for being master DOFS
DOFm_13  = [Indexes{1}.MASTER(:);Indexes{3}.MASTER(:)] ; 
DOFm_24  = [Indexes{2}.MASTER(:);Indexes{4}.MASTER(:)] ; 
DOFs_13  = [Indexes{1}.SLAVES(:);Indexes{3}.SLAVES(:)] ; 
DOFs_24  = [Indexes{2}.SLAVES(:);Indexes{4}.SLAVES(:)] ; 
 
ndofm_13 = length(DOFm_13 )  ;
ndofm_24 = length(DOFm_24)  ;
ndofm_cand = ndofm_13+ndofm_24 ; 
if ndofm_cand == size(Dcomp,2)
    DOFm = sort([DOFm_13;DOFm_24]) ; 
    DOFs = sort([DOFs_13; DOFs_24]  ) ; 
elseif ndofm_cand > size(Dcomp,2)
    error('Option not implemented (yet)')
else
    error('The number of candidats for being MASTER DOFs is less than the total number of domain modes')
end

Dcomp_m = Dcomp(DOFm,:) ; 

if rank(Dcomp_m)~=size(Dcomp,2)
    error('Choice of master DOFs conducive to rank-deficient Dcomp matrix')
end
