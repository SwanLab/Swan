classdef RHSIntegratorElasticMicro < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        material
        globalConnec

        fun
    end
    
    methods (Access = public)

        function obj = RHSIntegratorElasticMicro(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj,strainBase,test)
      %      oX     = zeros(obj.dim.ndimf,1);
     %       nVoigt = size(obj.material.evaluate(oX),1);
     %       basis   = diag(ones(nVoigt,1));
     %       Fvol = zeros(obj.dim.ndofs, nVoigt);
     %       for iVoigt = 1:nVoigt
      %          vstrain = basis(iVoigt,:);
                FvolE = obj.computeStrainRHS(strainBase);
       %         Fvol(:,iVoigt)  = obj.assembleVector(FvolE);
            Fvol = obj.assembleVector(FvolE,test);
       %     end
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

        function b = assembleVector(obj, forces,test)
        %    s.dim          = obj.dim;
        %    s.globalConnec = obj.globalConnec;
        %    s.nnodeEl      = []; % size(obj.geometry.dNdx,2);
%             F(:,1,:) = forces;
            assembler = AssemblerFun();
            b = assembler.assembleV(forces, test);
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
        
        
        function rhs = computeStrainRHS(obj,f)
            % quad    = Quadrature.create(obj.mesh, 1);
            % xV      = quad.posgp;
            % ngaus   = size(xV,2);
            % SymN = Grad(test);
            % symN  = SymN.evaluate(xV);            
            % dV      = obj.mesh.computeDvolume(quad);
            % C       = obj.material.evaluate(xV);
            % vstrain = strain.evaluate(xV);
            % nnodeE  = obj.mesh.nnodeElem;
            % ndim    = obj.mesh.ndim;
            % ndofE   = nnodeE*ndim;
            % nElem   = obj.mesh.nelem;
            % F       = zeros(ndofE,ngaus,nElem);
            % for i = 1:ndofE
            %     symTest  = symN(:,:,:,:);
            %     sig      = pagetensorprod(C,vstrain,[3 4],[1 2],4,2);
            %     df       = pagetensorprod(symTest,sig,[1 2],[1 2],2,2);
            %     F(i,:,:) = -df.*dV;
            % end

            quad    = Quadrature.create(obj.mesh, 2);

            nElem  = obj.mesh.nelem;
            rhs    = zeros(size(f,2),nElem);
            J = Jacobian(obj.mesh);
            detJ = Det(J);

            xV = quad.posgp;
            w  = quad.weigp;
            for i = 1:size(f,2)
                int = squeezeParticular((f{i}.*detJ)*w',2);
                rhs(i,:) = rhs(i,:) + int.evaluate(xV);
            end



        end

    end
    
end