classdef LHSintegrator_Mass_RT < LHSintegrator

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
            
            JGlob = obj.mesh.geometry.jacobian;
            Jdet = obj.mesh.geometry.detJ;
            
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
                                Jd = Jdet(:,igauss);
                                Ni = pagemtimes(squeeze(shapes(inode,igauss,:))',JGlob);
                                Nj = pagemtimes(squeeze(shapes(jnode,igauss,:))',JGlob);
                                v = squeeze(pagemtimes(Ni,pagetranspose(Nj)));
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*(dvol./Jd).*(sides(:,idof).*sides(:,jdof));
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