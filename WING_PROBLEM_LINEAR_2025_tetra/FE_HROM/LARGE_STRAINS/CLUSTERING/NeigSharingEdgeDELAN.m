function NeighboringElements = NeigSharingEdgeDELAN(CN,InvCNmatrix)

if nargin == 0
    load('tmp2.mat')
end

NeighboringElements = zeros(size(CN));
NumberNElements = zeros(size(CN,1),1) ;

    innodeGLO = [1:size(CN,2),1] ;


for ielem = 1:size(CN,1)
    
    inodeLOC = 1 ;
    while NumberNElements(ielem)  <=size(CN,2) && inodeLOC <=size(CN,2)
        innode = innodeGLO(inodeLOC) ;
        IndCommonElements = find(InvCNmatrix(CN(ielem,innode),:)>0 );
        CommonElements_i  =  InvCNmatrix(CN(ielem,innode),IndCommonElements) ;
        jnnode =  innodeGLO(inodeLOC+1) ; 
        IndCommonElements = find(InvCNmatrix(CN(ielem,jnnode),:)>0 );
        CommonElements_j  =  InvCNmatrix(CN(ielem,jnnode),IndCommonElements)  ;
        neig_elem = intersect(CommonElements_i,CommonElements_j) ;
        neig_elem = setdiff(neig_elem,ielem) ; 
        
        if ~isempty(neig_elem)
            III =  intersect(NeighboringElements(ielem,1:NumberNElements(ielem) ),neig_elem) ;
            if isempty(III)
                NumberNElements(ielem) = NumberNElements(ielem) + 1;
                NeighboringElements(ielem,NumberNElements(ielem)) = neig_elem ;
                
            end
        end
          inodeLOC = inodeLOC + 1;
        
    end
    
end
