classdef Anodal2gausComputer < handle

    properties (Access = public)
        A_nodal_2_gauss
    end

    properties (Access = private)
        nnode
        nelem
        npnod
        ngaus
        connec
        shape
    end

    methods (Access = public)
        function obj = Anodal2gausComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeA();
        end

        function intX = integrateP1FunctionWithShapeFunction(obj,cParams)
            ndof = size(obj.A_nodal_2_gauss{1},2);
            intX = zeros(ndof,1);
            for igaus = 1:obj.ngaus
                dVG  = cParams.dV(:,igaus);
                xG   = cParams.x(:,igaus);
                A    = obj.A_nodal_2_gauss{igaus};
                intX = intX + A'*(xG.*dVG);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.nnode  = cParams.nnode;
            obj.nelem  = cParams.nelem;
            obj.npnod  = cParams.npnod;
            obj.ngaus  = cParams.ngaus;
            obj.connec = cParams.connec;
            obj.shape  = cParams.shape;
        end

        function computeA(obj)
            A0    = sparse(obj.nelem,obj.npnod);
            A2g   = cell(obj.ngaus,1);
            nodes = obj.connec;
            for igaus = 1:obj.ngaus
                A2g{igaus} = A0;
                for inode = 1:obj.nnode
                    node   = nodes(:,inode);
                    shapeN = obj.shape(inode,igaus);
                    Ni = ones(obj.nelem,1)*shapeN;
                    A  = sparse(1:obj.nelem,node,Ni,obj.nelem,obj.npnod);
                    A2g{igaus} = A2g{igaus} + A;
                end
            end
            obj.A_nodal_2_gauss = A2g;
        end
    end

end