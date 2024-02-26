classdef Hexahedra_Linear < Interpolation

    methods (Access = public)

        function obj = Hexahedra_Linear(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function shape = computeShapeFunctions(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            shape = zeros(obj.nnode,ngaus,nelem);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            u = posgp(3,:,:);
            lcord = [-1 -1 -1;
                      1 -1 -1;
                      1  1 -1;
                     -1  1 -1;
                     -1 -1  1;
                      1 -1  1;
                      1  1  1;
                     -1  1  1];
            for inode = 1:obj.nnode
                shape(inode,:,:)=(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,2)*t).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
            end
        end

        function deriv = computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
            nelem = size(posgp,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            s = posgp(1,:,:);
            t = posgp(2,:,:);
            u = posgp(3,:,:);
            lcord = [-1 -1 -1;
                      1 -1 -1;
                      1  1 -1;
                     -1  1 -1;
                     -1 -1  1;
                      1 -1  1;
                      1  1  1;
                     -1  1  1];
            for inode = 1:obj.nnode
                deriv(1,inode,:,:) = lcord(inode,1).*(ones(1,size(posgp,2))+lcord(inode,2)*t).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
                deriv(2,inode,:,:) = lcord(inode,2).*(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
                deriv(3,inode,:,:) = lcord(inode,3).*(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,2)*t)/8;
            end
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.ndime = 3;
            obj.nnode = 8;
            % obj.isoDv = 8;
            obj.pos_nodes=[-1 -1 -1;
                +1 -1 -1;
                +1 +1 -1;
                -1 +1 -1;
                -1 -1 +1;
                +1 -1 +1;
                +1 +1 +1;
                -1 +1 +1];
        end

    end

end
