classdef LHSintegrator_Mass < LHSintegrator

    properties (Access = public)
        elemMass
    end
    
    properties (Access = private)
        field
    end

    
    methods (Access = public)

        function obj = LHSintegrator_Mass(cParams)
            %             obj.init(cParams);
            obj.mesh  = cParams.mesh;
            obj.field = cParams.field;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrixField(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            f = obj.field;
            shapes = f.interpolation.shape;
            quad   = f.quadrature;
            dvolu  = obj.mesh.computeDvolume(quad);
            ngaus  = f.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            %             nnode  = obj.mesh.nnodeElem;
            ndimf  = f.dim.ndimf;
            nnode  = f.dim.nnodeElem;

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
            M = zeros(nnode*ndimf,nnode*ndimf,nelem);
            dvolu = dvolu';
            for igauss = 1 :ngaus
                for inode= 1:nnode
                    for jnode= 1:nnode
                        for iunkn= 1:ndimf
                            for junkn= 1:ndimf
                                idof = ndimf*(inode-1)+iunkn;
                                jdof = ndimf*(jnode-1)+junkn;
                                dvol = dvolu(:,igauss);
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
            obj.elemMass = lhs;  
        end

        function lhs = assembleMatrixField(obj, Ae)
            s.dim          = obj.field.dim;
            s.globalConnec = obj.field.connec;
            s.nnodeEl      = obj.field.dim.nnodeElem;
            assembler = Assembler(s);
            lhs = assembler.assembleFields(Ae, obj.field, obj.field);
        end

    end

end