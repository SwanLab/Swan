classdef LHSintegrator_MassFun < handle

    properties (Access = private)
        mesh
        fun
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegrator_MassFun(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun      = cParams.fun;
            obj.mesh     = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.fun.order;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end
    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            shapes = obj.fun.computeShapeFunctions(obj.quadrature);
            dVolu  = obj.mesh.computeDvolume(obj.quadrature);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nDimf  = obj.fun.ndimf;
%             nDofs  = numel(obj.fun.fValues);
            nNodE  = obj.mesh.nnodeElem;
            nDofE  = nNodE*nDimf;

            % One dimension
            %             lhs = zeros(nnode,nnode,nelem);
            %             for igaus = 1:ngaus
            %                 dv(1,1,:) = dvolu(igaus,:);
            %                 Ni = shapes(:,igaus);
            %                 Nj = shapes(:,igaus);
            %                 NiNj = Ni*Nj';
            %                 Aij = bsxfun(@times,NiNj,dv);
            %                 lhs = lhs + Aij;
            %             end

            % N dimensions, pending optimization
            M = zeros(nDofE,nDofE,nElem);
            dVolu = dVolu';
            for igauss = 1 :nGaus
                for inode= 1:nNodE
                    for jnode= 1:nNodE
                        for iunkn= 1:nDimf
                            for junkn= 1:nDimf
                                idof = nDimf*(inode-1)+iunkn;
                                jdof = nDimf*(jnode-1)+junkn;
                                dvol = dVolu(:,igauss);
                                Ni = shapes(inode,igauss,:);
                                Nj = shapes(jnode,igauss,:);
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*dvol;
                            end
                        end
                    end
                end
            end
            lhs = M;

        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = obj.fun; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs);
        end

    end

end