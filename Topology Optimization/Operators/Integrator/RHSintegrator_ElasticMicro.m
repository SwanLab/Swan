classdef RHSintegrator_ElasticMicro < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        vstrain
        material
        geometry
        dvolume
        globalConnec
        quadrature
    end
    
    methods (Access = public)

        function obj = RHSintegrator_ElasticMicro(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createGeometry();
            obj.computeDvolume();
        end

        function Fext = compute(obj)
            Fvol       = obj.computeStrainRHS(obj.vstrain);
            forces     = squeeze(Fvol);
            FextSupVol = obj.assembleVector(forces);
            Fpoint     = obj.computePunctualFext();
            Fext = FextSupVol +  Fpoint;
        end

        function R = computeReactions(obj, K)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet;
            dirichV = bc.dirichlet_values;
            if ~isempty(dirich)
                R = -K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndofs(:)),1);
            end

        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim                = cParams.dim;
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.BC;
            obj.material           = cParams.material;
            obj.globalConnec       = cParams.globalConnec;
            obj.vstrain            = cParams.vstrain;
        end
       
        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
        end

        function createGeometry(obj)
            q = obj.quadrature;
            int = obj.mesh.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function computeDvolume(obj)
            q = obj.quadrature;
            obj.dvolume = obj.mesh.computeDvolume(q)';
        end

        function b = assembleVector(obj, forces)
            s.dim          = obj.dim;
            s.globalConnec = obj.globalConnec;
            s.nnodeEl      = []; % size(obj.geometry.dNdx,2);
%             F(:,1,:) = forces;
            assembler = Assembler(s);
            b = assembler.assembleV(forces);
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            neumann       = obj.boundaryConditions.neumann;
            neumannValues = obj.boundaryConditions.neumann_values;
            Fp = zeros(obj.dim.ndofs,1);
            if ~isempty(neumann)
                Fp(neumann) = neumannValues;
            end
        end
        
        function F = computeStrainRHS(obj,vstrain)
            Cmat  = obj.material.C;
            nunkn = obj.dim.ndimf;
            nstre = size(Cmat,1);
            nelem = size(Cmat,3);
            nnode = obj.dim.nnodeElem;
            ngaus = obj.quadrature.ngaus;

            eforce = zeros(nunkn*nnode,ngaus,nelem);
            sigma = zeros(nstre,ngaus,nelem);
            s.dim = obj.dim;
            s.geometry = obj.geometry;
            Bcomp = BMatrixComputer(s);
            for igaus = 1:ngaus
                Bmat    = Bcomp.compute(igaus);
                dV(:,1) = obj.dvolume(:,igaus);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij = squeeze(Cmat(istre,jstre,:,igaus));
                        vj  = vstrain(jstre);
                        si  = squeeze(sigma(istre,igaus,:));
                        sigma(istre,igaus,:) = si + Cij*vj;
                    end
                end
                for iv = 1:nnode*nunkn
                    for istre = 1:nstre
                        Biv_i = squeeze(Bmat(istre,iv,:));
                        si    = squeeze(sigma(istre,igaus,:));
                        Fiv   = squeeze(eforce(iv,igaus,:));
                        eforce(iv,igaus,:) = Fiv + Biv_i.*si.*dV;
                    end
                end
            end
            F = -eforce;
        end

    end
    
end