classdef Hexahedra_Linear < Interpolation

    methods (Access = public)

        function obj = Hexahedra_Linear(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

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

        function shape = evaluateShapeFunctions(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            shape = zeros(obj.nnode,ngaus,nelem);
            s = xV(1,:,:);
            t = xV(2,:,:);
            u = xV(3,:,:);
            lcord = [-1 -1 -1;
                      1 -1 -1;
                      1  1 -1;
                     -1  1 -1;
                     -1 -1  1;
                      1 -1  1;
                      1  1  1;
                     -1  1  1];
            for inode = 1:obj.nnode
                shape(inode,:,:)=(ones(1,size(xV,2))+lcord(inode,1)*s).*(ones(1,size(xV,2))+lcord(inode,2)*t).*(ones(1,size(xV,2))+lcord(inode,3)*u)/8;
            end
        end

        function deriv = evaluateShapeDerivatives(obj,xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            s = xV(1,:,:);
            t = xV(2,:,:);
            u = xV(3,:,:);
            lcord = [-1 -1 -1;
                      1 -1 -1;
                      1  1 -1;
                     -1  1 -1;
                     -1 -1  1;
                      1 -1  1;
                      1  1  1;
                     -1  1  1];
            for inode = 1:obj.nnode
                deriv(1,inode,:,:) = lcord(inode,1).*(ones(1,size(xV,2))+lcord(inode,2)*t).*(ones(1,size(xV,2))+lcord(inode,3)*u)/8;
                deriv(2,inode,:,:) = lcord(inode,2).*(ones(1,size(xV,2))+lcord(inode,1)*s).*(ones(1,size(xV,2))+lcord(inode,3)*u)/8;
                deriv(3,inode,:,:) = lcord(inode,3).*(ones(1,size(xV,2))+lcord(inode,1)*s).*(ones(1,size(xV,2))+lcord(inode,2)*t)/8;
            end
        end

    end

end
