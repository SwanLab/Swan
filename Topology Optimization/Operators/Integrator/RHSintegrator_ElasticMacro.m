classdef RHSintegrator_ElasticMacro < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        newBCs
    end
    
    methods (Access = public)

        function obj = RHSintegrator_ElasticMacro(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj)
%             Fsup       = obj.computeSuperficialFext;
%             Fvol       = obj.computeVolumetricFext;
%             forces     = squeeze(Fsup + Fvol);
%             forces     = squeeze(Fvol);
%             FextSupVol = obj.assembleVector(forces);
            Fext     = obj.computePunctualFext();
%             Fext = FextSupVol +  Fpoint;
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
            % obj.newBCs = cParams.newBCs;
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

    end
    
end