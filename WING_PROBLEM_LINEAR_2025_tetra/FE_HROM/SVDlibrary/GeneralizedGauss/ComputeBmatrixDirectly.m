function   [N,B] = ComputeBmatrixDirectly(COORnodes,xLOC,TypeElement)

if nargin == 0
    load('tmp1.mat')
end

ndim =  size(COORnodes,2) ;
nnode = size(COORnodes,1) ;
B = zeros(ndim,nnode) ;
N  = zeros(1,nnode)  ;


switch  TypeElement
    case 'Quadrilateral'
        if nnode == 9
            m = [2,2] ;
        elseif nnode == 4
            m = [1,1] ;
        else
            error('element not implemented')
        end
        
        
    case 'Hexahedra'         
        if nnode == 8
            m =[1,1,1] ;
        elseif nnode == 27
            m = [2,2,2] ;
        else
            error('Element not implemented ')
        end
        
        
    otherwise
        error('Element not implemented ')
end

SCALE_ACTIVE = 1; 
LelemSCALING = norm(COORnodes(2,:)-COORnodes(1,:)) ;
coorREF = COORnodes(1,:) ; 
if SCALE_ACTIVE == 1 
    % To avoid bad-scaling of the shape functions (coefficient matrix interpolation)
       COORnodes = bsxfun(@minus,COORnodes',coorREF')'/LelemSCALING   ;
       xLOC = (xLOC-coorREF)/LelemSCALING ; 
        factorDERIVATIVE = 1/LelemSCALING ; 
else
    factorDERIVATIVE = 1;
end


[Nloc,Bloc]  =  PolynomialInterpolationShapeFunctions(COORnodes,m ,xLOC) ;
N = Nloc ;

for idim = 1:length(Bloc)
    B(idim,:) = Bloc{idim}*factorDERIVATIVE ;
end


