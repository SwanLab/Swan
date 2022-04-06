classdef LHSintegrator_Stiffness < LHSintegrator
    
    properties (Access = private)
        geometry
    end

    methods (Access = public)
        
        function obj = LHSintegrator_Stiffness(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
        
    end
    
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
            dvolu = obj.mesh.computeDvolume(obj.quadrature);
            ngaus = obj.quadrature.ngaus;
            nelem = obj.mesh.nelem;
            nstre = obj.dim.nstre;
            ndpe  = obj.dim.ndofPerElement;
            lhs = zeros(ndpe,ndpe,nelem);
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
                for istre = 1:nstre
                    BmatI = Bmat(istre,:,:);
                    BmatJ = permute(BmatI,[2 1 3]);
                    dNdN = bsxfun(@times,BmatJ,BmatI);
                    dv(1,1,:) = dvolu(igaus, :);
                    inc = bsxfun(@times,dv,dNdN);
                    lhs = lhs + inc;
                end
            end
        end
        
   end
    
   methods (Access = private)
       
        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function Bcomp = createBComputer(obj)
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = obj.globalConnec;
%             s.dofsInElem   = obj.dofsInElem;
            Bcomp = BMatrixComputer(s);
        end
       
   end
    
end