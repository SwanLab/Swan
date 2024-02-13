classdef RHSintegrator_ElasticMicro < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        vstrain
        material
        dvolume
        globalConnec
        quadrature

        fun
    end
    
    methods (Access = public)

        function obj = RHSintegrator_ElasticMicro(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeDvolume();
        end

        function Fext = compute(obj)
            nVoigt = size(obj.material.evaluate([0;0]),1);
            basis   = diag(ones(nVoigt,1));
            Fvol = zeros(obj.dim.ndofs, nVoigt);
            for iVoigt = 1:nVoigt
                vstrain = basis(iVoigt,:);
                FvolE = obj.computeStrainRHS(vstrain);
                Fvol(:,iVoigt)  = obj.assembleVector(FvolE);
            end
            Fpoint = obj.computePunctualFext();
            Fext = Fvol + Fpoint;
        end

        function R = computeReactions(obj, K)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet_dofs;
            dirichV = bc.dirichlet_vals;
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
            obj.fun                = cParams.fun;
        end
       
        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            obj.quadrature = q;
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
            s.fun = obj.fun;
            assembler = AssemblerFun(s);
            b = assembler.assembleV(forces, obj.fun);
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            neumann       = obj.boundaryConditions.pointload_dofs;
            neumannValues = obj.boundaryConditions.pointload_vals;
            Fp = zeros(obj.dim.ndofs,1);
            if ~isempty(neumann)
                Fp(neumann) = neumannValues;
            end
        end
        
        
        function F = computeStrainRHS(obj,vstrain)
            xV    = obj.quadrature.posgp;
            Cmat  = obj.material.evaluate(xV);
            nunkn = obj.dim.ndimf;
            nstre = size(Cmat,1);
            nelem = size(Cmat,3);
            nnode = obj.dim.nnodeElem;
            ngaus = obj.quadrature.ngaus;

            eforce = zeros(nunkn*nnode,ngaus,nelem);
            sigma = zeros(nstre,ngaus,nelem);

            a.mesh       = obj.mesh;
            a.fValues    = sigma;
            a.quadrature = obj.quadrature;
            sigmaF = FGaussDiscontinuousFunction(a);

            sigmaF.ndimf = size(obj.mesh.coord,2);
            s.fun  = sigmaF;
            s.dNdx = sigmaF.computeCartesianDerivatives(obj.quadrature);

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