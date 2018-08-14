function [ nproc,coeff ] = nprocedure( etype,nnode )

% Remark : lumped matrix is computed according to the comments described in
% FEM its basis & fund,(Zienkiewicz,taylos,zhu), pp 568, vol I 
% nproc = 1,2,3; Row sum procedure, Diagonal scaling, Quadrature using nodal points
coeff = [];
switch etype
    case 'LINEAR_TRIANGLE'
        switch nnode
            case 3
                nproc = 1;
        end
    case {'QUAD'}
        switch nnode
            case 4
                nproc = 1;
            case 8
                nproc = 2;
                a=1/36; b=8/36; 
                coeff(1)=a;
                coeff(2)=a;
                coeff(3)=a;
                coeff(4)=a;
                coeff(5)=b;
                coeff(6)=b;
                coeff(7)=b;
                coeff(8)=b;
            case 9
                nproc = 1;
        end
end

end

