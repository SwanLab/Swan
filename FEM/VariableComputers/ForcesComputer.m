classdef ForcesComputer < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
        vstrain
        material
        geometry
        dvolume
        dofsInElem
        globalConnec
    end
    
    methods (Access = public)

        function obj = ForcesComputer(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj)
            Fsup       = obj.computeSuperficialFext;
            Fvol       = obj.computeVolumetricFext;
            forces     = squeeze(Fsup + Fvol);
            FextSupVol = obj.assembleVector(forces);
            Fpoint     = obj.computePunctualFext();
            Fext = FextSupVol +  Fpoint;
        end

        function R = computeReactions(obj, K)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet{1};
            dirichV = bc.dirichlet_values{1};
            if ~isempty(dirich)
                R = -K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndof(:)),1);
            end

        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim                = cParams.dim;
            obj.mesh               = cParams.mesh;
            obj.boundaryConditions = cParams.BC;
            obj.material           = cParams.material;
            obj.geometry           = cParams.geometry;
            obj.dvolume            = cParams.dvolume';
            obj.globalConnec       = cParams.globalConnec;
            if isfield(cParams, 'vstrain')
                obj.vstrain = cParams.vstrain;
            end
        end

        function Fs = computeSuperficialFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            ndimf = d.ndimField;
            nelem = d.nelem;
            Fs = zeros(nnode*ndimf,1,nelem);
        end
        
        function Fv = computeVolumetricFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            ndimf = d.ndimField;
            nelem = d.nelem;
            Fv = zeros(nnode*ndimf,1,nelem);
            if ~isempty(obj.vstrain)
                Fv = Fv + obj.computeStrainRHS(obj.vstrain);
            end
        end

        function b = assembleVector(obj, forces)
            s.dim          = obj.dim;
            s.globalConnec = obj.globalConnec;
            assembler = Assembler(s);
            b = assembler.assembleV(forces);
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            neumann       = obj.boundaryConditions.neumann;
            neumannValues = obj.boundaryConditions.neumann_values;
            Fp = zeros(obj.dim.ndof,1);
            if ~isempty(neumann)
                Fp(neumann) = neumannValues;
            end
        end

        
        function F = computeStrainRHS(obj,vstrain)
%             Cmat  = obj.material.C(:,:,1);
            Cmat  = obj.material.C;
            ngaus = obj.dim.ngaus;
            nunkn = obj.dim.ndimField;
            nstre = obj.dim.nstre;
            nelem = obj.dim.nelem;
            nnode = obj.dim.nnode;

            eforce = zeros(nunkn*nnode,ngaus,nelem);
            sigma = zeros(nstre,ngaus,nelem);
            s.dim = obj.dim;
            s.geometry = obj.geometry;
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
            for igaus = 1:ngaus
                Bmat    = Bcomp.computeBmat(igaus);
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