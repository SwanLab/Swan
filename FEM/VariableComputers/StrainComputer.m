classdef StrainComputer < handle
    
    properties (Access = private)
        dim
        mesh
        connec
        geometry
        quadrature
        displacement
        dispField
    end
    
    methods (Access = public)
        
        function obj = StrainComputer(cParams)
            obj.init(cParams)
%             obj.createGeometry();
        end

        function strain = compute(obj)
            strain = obj.computeStrain();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dispField          = cParams.dispField;
            obj.dim                = cParams.dispField.dim;
            obj.mesh               = cParams.mesh;
            obj.quadrature         = cParams.quadrature;
            obj.displacement       = cParams.displacement;
            obj.connec             = cParams.dispField.connec;
            obj.geometry           = cParams.dispField.geometry;
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
        
        function strain = computeStrain(obj)
            d_u = obj.displacement;
            connec = obj.connec;
            nelem   = size(connec,1);
            nstre   = obj.getNstre();
            nnodeEl = obj.dim.nnodeElem;
            nunkn   = obj.dim.ndimf;
            ngaus   = obj.quadrature.ngaus;
            strain = zeros(nstre,nelem,ngaus);
            for igaus = 1:ngaus
                Bmat = obj.computeB(igaus);
                for istre=1:nstre
                    for inode=1:nnodeEl
                        nodes = connec(:,inode);
                        for idime = 1:nunkn
                            dofs = nunkn*(nodes - 1) + idime;
                            ievab = nunkn*(inode-1)+idime;
                            B = squeeze(Bmat(istre,ievab,:));
                            u = d_u(dofs);
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(B.*u)';
                        end
                        
                    end
                end
            end
            strain = permute(strain, [3 1 2]);
        end

        function nstre = getNstre(obj)
            nstreVals = [2, 3, 6];
            nstre = nstreVals(obj.dim.ndimf);
        end

        function Bmat = computeB(obj,igaus)
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
            Bmat = Bcomp.computeBmat(igaus);
        end
    end
    
end