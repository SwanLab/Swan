function    invJt= inverseTRANSvectorize(Jinp,ndim,detJinp)
% Given a matrix Jinp = [Jinp_1; Jinp_2 ... ; Jinp_nelem]  (Jinp_i is ndim x ndim),
% inverseTRANSvectorize returns a matrix consisting of the  transponse of
% the inverse of each block matrix
%  invJt = [inv(Jinp_1'); inv(Jinp_2') ... ;inv(Jinp_nelem') ]
% ndim = 2,3
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
if nargin == 0
    J1 = [3 5 6; 7 8 7 ; 0 0 3] ;  d1 =inv(J1') ;
    J2 = 4*[3 5 6; 7 8 7 ; 0 0 3] ;  d2 =inv(J2') ;
    J3 = 6*[3 5 6; 7 8 7 ; 0 0 3] ;  d3 =inv(J3') ;
    Jinp = [J1; J2; J3]; invJt_check = [d1;d2;d3] ;
    detJinp = [] ;
    ndim = 3 ;
end

if isempty(detJinp)
     detJinp= determinantVECTORIZE(Jinp,ndim) ; 
end

%
J = cell(ndim,ndim) ;
for i=1:ndim
    iglo = i:ndim:size(Jinp,1) ;
    for j=1:ndim
        jglo = j ;
        J{i,j} = Jinp(iglo,jglo) ;
    end
end
% ------------------------------------------------
    invJt = zeros(size(detJinp)) ;

if ndim ==2
      % i = 1
    % -----
    i= 1; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJt(iglo,1) = J{2,2}./detJinp ;
    invJt(iglo,2) = - J{2,1}./detJinp ;
    % i = 2
    % -----
    i= 2; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJt(iglo,1) =  -J{1,2}./detJinp ;
    invJt(iglo,2) =   J{1,1}./detJinp ;
    
elseif ndim == 3
    % i = 1
    % -----
    i= 1; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJt(iglo,1) = (J{2,2}.*J{3,3} - J{2,3}.*J{3,2})./detJinp ;
    invJt(iglo,2) = (J{2,3}.*J{3,1} - J{2,1}.*J{3,3})./detJinp ;
    invJt(iglo,3) = (J{2,1}.*J{3,2} - J{2,2}.*J{3,1})./detJinp ;
    % i = 2
    % -----
    i= 2; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJt(iglo,1) = (J{1,3}.*J{3,2} - J{1,2}.*J{3,3} )./detJinp ;
    invJt(iglo,2) = (J{1,1}.*J{3,3} - J{1,3}.*J{3,1} )./detJinp ;
    invJt(iglo,3) = (J{1,2}.*J{3,1} - J{1,1}.*J{3,2} )./detJinp ;
    % i = 3
    % -----
    i= 3; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJt(iglo,1) = (J{1,2}.*J{2,3} - J{1,3}.*J{2,2} )./detJinp ;
    invJt(iglo,2) = ( J{1,3}.*J{2,1} - J{1,1}.*J{2,3} )./detJinp ;
    invJt(iglo,3) = (J{1,1}.*J{2,2} - J{1,2}.*J{2,1})./detJinp ;
end


