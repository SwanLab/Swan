function [N,BeTILDE,coor_PARENTDOMAIN] = ShapeFunDer_inversemapping(Xe,xLOC,TypeElement,DATAinp)

if nargin == 0
    load('tmp.mat')
    DATAinp.ORDER_INVERSE_ELEMENT_inverse_mapping = 2; %
    DATAinp.NPOINTS_ONE_DIRECTION_inverse_mapping  =3 ;
     DATAinp.InterpolationMethod  = 'FEnodes'; %(DATA,'InterpolationMethod','ScatteredDataMatlab') ;

end
% 3-may-2020. Inverse Mapping method
% Xe = COOR(NODESloc,:)' ;
nnodeE = size(Xe,2)  ;
DATAinp = DefaultField(DATAinp,'ipoint',1) ;
switch TypeElement
    case 'Hexahedra'
        if nnodeE== 27
            
            disp(['Performing inverse mapping ...',' ielem = ',num2str(DATAinp.ipoint)]) ;
            xiNEW =  InverseMappingHexag27(Xe',xLOC,DATAinp) ;
        elseif nnodeE== 20
            xiNEW =  InverseMappingHexag20(Xe',xLOC,DATAinp) ;
        else
            error('Not implemented...')
        end
    case 'Quadrilateral'
        if nnodeE== 9
            
            %disp(['Performing inverse mapping ...',' ielem = ',num2str(DATAinp.ipoint)]) ;
            xiNEW =  InverseMappingQuad9(Xe',xLOC,DATAinp) ;
        else
            nnodeE
            error('Not implemented...')
        end
        
    otherwise
        error('Not implemented...')
        
end



switch TypeElement
    case 'Hexahedra'
        
        if   nnodeE== 20
            
            COORparent  = COOR_Hexag20 ;
            ORDER_POLYNOMIALS=[2,2,2] ;
            DATAshape = ShapeFunCoefficientsSEREN(COORparent,ORDER_POLYNOMIALS) ;
            DATAlocSHAPE.DATAshape  = DATAshape;
            xLIM = [] ;
            DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
            [N, ~,~ ]=    ShapeFunctionFEseren(xLIM,xiNEW',DATAlocSHAPE) ;
            
            BeTILDE = [] ;
            coor_PARENTDOMAIN = [] ;
        elseif nnodeE == 27
            coor_PARENTDOMAIN  = xiNEW ;
            [N,dershapef] = ComputeElementShapeFun_givenpoint(TypeElement,nnodeE,[],xiNEW) ;
            
            
            
            
            BeXi = dershapef(:,:,1) ;
            % Jacobian Matrix
            Je = Xe*BeXi' ;
            BeTILDE = inv(Je)'*BeXi ;
            
        else
            error('Option not implemented')
            
        end
        
    otherwise
        
        coor_PARENTDOMAIN  = xiNEW ;
        [N,dershapef] = ComputeElementShapeFun_givenpoint(TypeElement,nnodeE,[],xiNEW) ;
        BeXi = dershapef(:,:,1) ;
        % Jacobian Matrix
        Je = Xe*BeXi' ;
        BeTILDE = inv(Je)'*BeXi ;
end


disp(['Error inverse mapping'])
% New coordinates
Xe_MAPPED = N*Xe' ;
DIFF_MAPPING = norm(Xe_MAPPED-xLOC)/norm(Xe)*100 ;
disp(['Error mapping =',num2str(DIFF_MAPPING),' %'])
disp('')
