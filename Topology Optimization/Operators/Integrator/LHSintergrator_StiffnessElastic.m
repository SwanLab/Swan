classdef LHSintergrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        geometry
    end

    methods (Access = public)

        function obj = LHSintergrator_StiffnessElastic(cParams)
            obj.init(cParams);
            obj.material = cParams.material;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs   = obj.computeElementalLHS();
            LHS   = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS2(obj)
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
                    for jstre = 1:nstre
                        Cij   = squeeze(Cmat(istre,jstre,:));
                        c(1,1,:) = Cij.*dV;
                        Bi = Bmat(istre,:,:);
                        CB = bsxfun(@times,Bi,c);
                        Bj = permute(Bmat(jstre,:,:),[2 1 3]);
                        t = bsxfun(@times,CB,Bj);
                        lhs = lhs + t;
                    end
                end
            end
        end

        %% From Element_DiffReact
%         function Ke = computeElementalStiffnessMatrix(obj)
%             obj.quadrature.computeQuadrature('LINEAR');
%             obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
%             ndof  = obj.dof.nunkn*obj.nnode;
%             ngaus = obj.quadrature.ngaus;
%             Ke = zeros(ndof,ndof,obj.nelem);
%             for igaus = 1:ngaus
%                 dShapeDx = obj.geometry.dNdx(:,:,:,igaus);
%                 Bmat     = obj.computeB(obj.dof.nunkn,obj.nelem,obj.nnode,dShapeDx);
%                 for istre = 1:obj.nstre
%                     BmatI = Bmat(istre,:,:);
%                     BmatJ = permute(Bmat(istre,:,:),[2 1 3]);
%                     dNdN = bsxfun(@times,BmatJ,BmatI);
%                     dv(1,1,:) = obj.geometry.dvolu(:,igaus);
%                     inc = bsxfun(@times,dv,dNdN);
%                     Ke = Ke + inc;
%                 end
%             end
%         end

        %% New attempt


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
                    for jstre = 1:nstre
                        Cij   = squeeze(Cmat(istre,jstre,:));
                        c(1,1,:) = Cij.*dV;
                        Bi = Bmat(istre,:,:);
                        CB = bsxfun(@times,Bi,c);
                        Bj = permute(Bmat(jstre,:,:),[2 1 3]);
                        t = bsxfun(@times,CB,Bj);
                        lhs = lhs + t;
                    end
                end
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