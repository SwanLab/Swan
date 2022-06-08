classdef LHSintegrator_StiffnessElastic < LHSintegrator

    properties (Access = public)
        elemStiff
    end
    
    properties (Access = private)
        geometry
        material
        field
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessElastic(cParams)
            obj.mesh = cParams.mesh;
            obj.field = cParams.field;
            obj.material = cParams.material;
          %  obj.interpolation = cParams.interpolation;
             obj.createQuadrature();
             obj.createInterpolation();
            obj.quadrature.computeQuadrature(obj.interpolation.order)
            obj.createGeometry();
        end

        function LHS = compute(obj)
%             disp('Elemental')
%             tic
%                 lhs = obj.computeElementalLHS();
%             toc
%             disp('Pagemtimes')
%             tic
              lhs   = obj.computeElementalLHSPagemtimes();
              obj.elemStiff = lhs;
%             toc
%             LHS = obj.assembleMatrix(lhs);
            LHS = obj.assembleMatrixField(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            C = obj.material.C;
            dvolu  = obj.mesh.computeDvolume(obj.field.quadrature);
            nstre  = size(C,1);
            nelem  = size(C,3);
            ngaus  = obj.field.quadrature.ngaus;
            npe    = obj.field.dim.ndofsElem;
            lhs = zeros(npe,npe,nelem);
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
                Cmat = C(:,:,:,igaus);
                dV    = dvolu(igaus,:)';
                for istre = 1:nstre
                    Bi = Bmat(istre,:,:);
                    for jstre = 1:nstre
                        Cij   = squeeze(Cmat(istre,jstre,:));
                        c(1,1,:) = Cij.*dV;
                        CB = bsxfun(@times,Bi,c);
                        Bj = permute(Bmat(jstre,:,:),[2 1 3]);
                        t = bsxfun(@times,CB,Bj);
                        lhs = lhs + t;
                    end
                end
            end
        end

        function lhs = computeElementalLHSPagemtimes(obj)
            dvolu  = obj.mesh.computeDvolume(obj.field.quadrature);
            ngaus  = obj.field.quadrature.ngaus;
            nelem  = size(obj.material.C,3);
            npe    = obj.dim.ndofsElem;
            % npe    = obj.dim.nnodeEl*obj.dim.ndimf;
            lhs = zeros(npe,npe,nelem);
            Bcomp = obj.createBComputer();
            Cmat = obj.material.C(:,:,:,1);
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
%                 Cmat = obj.material.C(:,:,:,igaus);
                dV(1,1,:) = dvolu(igaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, Bmat);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end

    end

    methods (Access = private)

        function Bcomp = createBComputer(obj)
            s.dim          = obj.field.dim;
            s.geometry     = obj.field.geometry;
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
        end

        function LHS = assembleMatrixField(obj, lhs)
            s.dim          = obj.field.dim;
            s.globalConnec = obj.field.connec;
            s.nnodeEl      = obj.field.dim.nnodeElem;
            assembler = Assembler(s);
            LHS = assembler.assembleFields(lhs, obj.field, obj.field);
        end
    end

end
