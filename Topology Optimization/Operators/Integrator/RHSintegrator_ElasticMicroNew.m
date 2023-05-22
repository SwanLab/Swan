classdef RHSintegrator_ElasticMicroNew < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        vstrain
        material
        dvolume
        globalConnec
        quadrature
    end
    
    methods (Access = public)

        function obj = RHSintegrator_ElasticMicroNew(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeDvolume();
        end

        function Fext = compute(obj)
            nVoigt = obj.material.nstre;
            basis   = diag(ones(nVoigt,1));
            Fvol = zeros(obj.dim.ndofs, nVoigt);
            for iVoigt = 1:nVoigt
                vstrain = basis(iVoigt,:);
                FvolE = squeeze(obj.computeStrainRHS(vstrain));
                Fvol(:,iVoigt)  = obj.assembleVector(FvolE);
            end
%             Fvol       = obj.computeStrainRHS(obj.vstrain);
%             forces     = squeeze(Fvol);
%             FextSupVol = obj.assembleVector(forces);
            Fpoint     = obj.computePunctualFext();
%             Fext = FextSupVol +  Fpoint;
            Fext = Fvol + Fpoint;
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
%             obj.vstrain            = cParams.vstrain;
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
            ngaus = obj.quadrature.ngaus;
            for igaus = 1:ngaus
                sigma  = obj.computeStress(vstrain);
                eforce = obj.computeEForce(sigma,igaus);
            end
            F = -eforce;
        end

        function sigmaFun = computeStress(obj, vstrain)
            Cmat  = obj.material.C;
            nElem = size(Cmat,3);
            vStr  = repmat(vstrain', [1  1 nElem]);
            sigma = pagemtimes(Cmat, vStr);

            a.mesh       = obj.mesh;
            a.fValues    = sigma;
            a.quadrature = obj.quadrature;
            sigmaFun = FGaussDiscontinuousFunction(a);
        end

        function eforce = computeEForce(obj, sigma, igaus)
            sigma.ndimf = size(obj.mesh.coord,2); 
            s.fun  = sigma;
            s.dNdx = sigma.computeCartesianDerivatives(obj.quadrature);
            Bcomp = BMatrixComputer(s);
            Bmat    = Bcomp.compute(igaus);
            
            dV(1,1,:) = obj.dvolume(:,igaus);
            Bok = permute(Bmat, [2 1 3]);
            Bsig = pagemtimes(Bok, sigma.fValues);
            eforce = pagemtimes(Bsig,dV);
        end

    end
    
end