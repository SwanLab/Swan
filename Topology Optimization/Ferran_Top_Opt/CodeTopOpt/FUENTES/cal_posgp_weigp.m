function [posgp,weigp,ngaus] = cal_posgp_weigp(type,ndime,nnode,ngaus)

switch type
    case {'TRIANGLE'} % SIMPLICIAL, TRIANG, 3 NODES
       switch nnode
            case 3
                %ngaus = 13; % 1 3 4 6 7 13
                ngaus = 1;
            case 6  
               
       end
       [posgp,weigp] = trian_gauss_const(ndime,ngaus);
    case 'QUAD' % quadrilateral
        switch nnode
            case 4
                ngaus = 4;
            case {8,9}
                ngaus = 9;
        end
        [posgp,weigp] = quad_gauss_const(ndime,ngaus);
    case 'HEXAHEDRA' % HEXAHEDRA 8 nodes
        ngaus = 8;
        [posgp,weigp] = quad_gauss_const(ndime,ngaus);
    otherwise
        error('undefined element type')
end

end

