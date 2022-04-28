classdef LHSintegrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        geometry
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessElastic(cParams)
            obj.init(cParams);
            obj.material = cParams.material;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
%             disp('Elemental')
%             tic
%                 lhs   = obj.computeElementalLHS();
%             toc
%             disp('Pagemtimes')
%             tic
                lhs   = obj.computeElementalLHSPagemtimes();
%             toc
            LHS   = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.dim.nelem;
            nstre  = obj.dim.nstre;
            npe    = obj.dim.ndofPerElement;
            lhs = zeros(npe,npe,nelem);
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
                Cmat = obj.material.C(:,:,:,igaus);
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
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.dim.nelem;
            npe    = obj.dim.ndofPerElement;
            lhs = zeros(npe,npe,nelem);
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
                Cmat = obj.material.C(:,:,:,igaus);
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
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = obj.globalConnec;
            Bcomp = BMatrixComputer(s);
        end

        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

    end

end