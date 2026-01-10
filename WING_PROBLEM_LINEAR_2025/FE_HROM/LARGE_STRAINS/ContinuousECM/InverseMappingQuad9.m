function xiNEW =  InverseMappingQuad9(COOR,xNEW,DATA)

if nargin == 0
    load('tmp1.mat')
    
end

DATA = DefaultField(DATA,'InterpolationMethod','ScatteredDataMatlab') ;
DATA = DefaultField(DATA,'ORDER_INVERSE_ELEMENT',2) ;
DATA = DefaultField(DATA,'NPOINTS_ONE_DIRECTION',40) ;

switch DATA.InterpolationMethod
    case 'FEnodes'
        nPOINTS = DATA.ORDER_INVERSE_ELEMENT +1;
    case 'ScatteredDataMatlab'
        
        nPOINTS = DATA.NPOINTS_ONE_DIRECTION ;
end

% NEW SET OF NODES IN THE PARENT DOMAIN
XX  = linspace(-1,1,nPOINTS) ; YY = XX;
[xx,yy ] = meshgrid(XX,YY) ;
xx = xx(:) ; yy = yy(:) ;
COORparent_rich =  [xx,yy];
%%  Physical coordinates corresponding to COORparent_rich
axy = 1;
COORparent  = COOR_Quad_n9(axy) ;
n=[2,2] ;
[PX]= CoordinateMatrixPolynomial(COORparent,n) ;
[Px]= CoordinateMatrixPolynomial(COORparent_rich,n) ;
COOR_rich = Px*(PX\COOR) ;

%%%  Now we create a higher-order element
switch DATA.InterpolationMethod
    case 'FEnodes'
        n = DATA.ORDER_INVERSE_ELEMENT*ones(1,2) ;
        [PX]= CoordinateMatrixPolynomial(COOR_rich,n) ;
        %
        [Px]= CoordinateMatrixPolynomial(xNEW,n) ;
        xiNEW = (Px*(PX\COORparent_rich))' ;
    case 'ScatteredDataMatlab'
        
        Finterp = scatteredInterpolant(COOR_rich,COORparent_rich(:,1),'natural') ;
        x = Finterp(xNEW) ;
        Finterp.Values = COORparent_rich(:,2);
        y = Finterp(xNEW) ;
        
        xiNEW  = [x,y]' ;
        %
        %
        %         for imode = 1:size(PHI,2)
        %     if imode == 1
        %         Finterp = scatteredInterpolant(COORg(:,1),COORg(:,2),PHI_OLD(:,imode),InterpolationMethod,ExtrapMethod) ;
        %     else
        %         Finterp.Values = PHI_OLD(:,imode) ;
        %     end
        %
        %     PHI(INDinterpolate,imode) = Finterp(COORsort(INDinterpolate,1),COORsort(INDinterpolate,2)) ;
        % end
        %
        
        
        
end
