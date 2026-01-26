function [IDXelem,IDXgauss ]= EquivalenceBetweenMeshes(COOR,IDX,CN,CNref,TypeIntegrand,...
    TypeElement)
% See Implementation.pdf, Collecting stress snapshots section.
% Equivalence between reference mesh ,CNref, and another mesh
%  (COOR,CN), where COOR = CNref(IDX,:)
if nargin == 0
    COOR = [0,2
        1,2
        2,2
        0,1
        1,1
        2,1
        0,0
        1,0
        2,0];

    
    CNref = [4 7 8 5
        2 5 6 3
        8 9 6 5
        1 4 5 2] ;
    IDX =[3 10 11 6  12 13 9 14 15]';
    CN = [  10 3  6 12
        10 12 13 11
        9 14 12 6
        15 13 12 14
        ];
    % ANSWER  IDXelem = [3 2 4 1]' ;
    TypeIntegrand = 'K';
    TypeElement = 'Quadrilateral' ;
    
    
end
[nnode,ndim] = size(COOR) ;
[nelem nnodeE ]= size(CN) ;
%-------------
IDXall = zeros(max(IDX),1) ;
IDXall(IDX) = 1:length(IDX) ;
CNtransf = IDXall(CN(:));
CNtransf = reshape(CNtransf,[],nnodeE) ;
% Now we have to compare row by row CNtransf with CNref. We recourse to the
% computation of the pseudo-center of gravity of each element
CGref = zeros(nelem,ndim) ;
for innodeE = 1:nnodeE
    CGref = CGref + COOR(CNref(:,innodeE),:) ;
end
CG = zeros(nelem,ndim) ;
for innodeE = 1:nnodeE
    CG = CG + COOR(CNtransf(:,innodeE),:) ;
end

IDXelem = knnsearch(CG,CGref);
% And what about connectivities per elemenent ? They may be different, so it is important
% to find also the mapping matrix, isn't it ? Or is it better to directly
% calculate the position of the Gauss points ? I think it is better to
% directly calculate the position of the Gauss points  and then use
% knnsearch element by element.
CNrefNS = CNtransf(IDXelem,:);  % Not-sorted reference matrix
[weig,posgp,N,dershapef] = ...
    ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
ngaus = length(weig) ;

% Loop over elements (it would be convenient to develop the vectorized counterpart...)
IDXnodLOC = zeros(size(CN)) ;
for ielem = 1:nelem
    NODESref = CNref(ielem,:) ;
    xREF = COOR(NODESref,:) ;
    NODES    = CNrefNS(ielem,:) ;
    x = COOR(NODES,:) ;
    % Gauss points position for xREF
    xGAUSSref = zeros(ngaus,ndim) ;
    xGAUSS = zeros(ngaus,ndim) ;
    for idim  =1:ndim
        xGAUSSref(:,idim) = N*xREF(:,idim) ;
        xGAUSS(:,idim) = N*x(:,idim) ;
    end
    
    IDXnodLOC(ielem,:) = knnsearch(xGAUSS,xGAUSSref)' ;
    
    
end

%% IDXgauss
% ---------
IDXgauss = zeros(nelem,ngaus) ;
for igaus = 1:ngaus
    IDXgauss(:,igaus) =  (IDXelem-1)*ngaus + IDXnodLOC(:,igaus) ; 
end

IDXgauss = reshape(IDXgauss',[],1) ; 



% for igaus = 1:ngaus
%     IDXgauss(igaus,:) = ngaus*(IDXelem'-1)+igaus ; 
% end
% IDXgauss = IDXgauss' ; 

 


end





%
%
%
% % Nonvectorized form
% IDXelemLOC = zeros(nelem,nnodeE) ;
% for ielem = 1:nelem
%     for innode = 1:nnodeE
%     IDXelemLOC(ielem,innode) = find(CNrefNS(ielem,innode)==CNref(ielem,:)) ;
%     end
% end
%
%
