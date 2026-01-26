function COORphysic = QuadrilateralFEshape_arbitraryORDER(OrderFEshapefunction,Xelem)



%We take as parent domain the biunit square [-1,1] x
% [-1,1], whose nodes are
xLOC = linspace(-1,1,OrderFEshapefunction+1) ; yLOC = xLOC ; 
[xx,yy]  = meshgrid(xLOC,yLOC) ;
xx = xx(:) ; yy = yy(:);
COORparent= [xx'; yy'] ;
% Now it proves convenient to order the columns of COORparent so that the
% first 4 columns are the coordinates of the vertices of the quadrilateral 
VERTICES = [-1 +1 +1 -1 ; -1 -1 +1 +1] ; % These are the coordinates of the 4 vertices 
mnonvertices = size(COORparent,2)-4;  
COORnon = zeros(size(COORparent,1),mnonvertices) ;  
inon = 0; 
for inode = 1:size(COORparent,2)
    COLUMNloc = COORparent(:,inode) ; 
    COMPARATION = bsxfun(@minus,VERTICES,COLUMNloc) ; 
    COMPARATION = sum(abs(COMPARATION),1) ; 
    ESCERO = find(COMPARATION==0) ; 
    if isempty(ESCERO)
        inon = inon + 1; 
        COORnon(:,inon) = COLUMNloc ; 
    end
end
COORparent = [VERTICES,COORnon]  ; 


% Now we wish to determine the positions of the points defined by the
% parent coordinates COORparent in the physical domain.  We do that by
% linear shape functions 
COORphysic = zeros(size(COORparent)); 
for innode = 1:size(COORparent,2)
    Ne = Quadrilateral4N(COORparent(:,innode)) ;
    for idim = 1:size(COORphysic,1)
        COORphysic(idim,innode) = Ne* Xelem(idim,:)' ; 
    end
end