classdef LHSintegrator_MassVect < LHSintegrator

    methods (Access = public)

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            shapes = obj.fun.computeShapeFunctions(obj.quadrature);
            dVolu  = obj.mesh.computeDvolume(obj.quadrature);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);
            nDimf  = obj.fun.ndimf;
            nNodE  = size(shapes,1);
            nDofE  = nNodE*nDimf;
            
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            sides = zeros(obj.mesh.nelem,obj.mesh.edges.nEdgeByElem);
            for ielem=1:obj.mesh.nelem
                sides(ielem,:) = ones(1,obj.mesh.edges.nEdgeByElem)-2.*(locPointEdge(ielem,:)~=1:obj.mesh.edges.nEdgeByElem);
            end
            
            J = [-1 -1;1 0];
            detJ = 1;

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
                                Ni = J*squeeze(shapes(inode,igauss,:));
                                Nj = J*squeeze(shapes(jnode,igauss,:));
                                v = squeeze(Ni'*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*dvol.*(sides(:,idof).*sides(:,jdof));
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