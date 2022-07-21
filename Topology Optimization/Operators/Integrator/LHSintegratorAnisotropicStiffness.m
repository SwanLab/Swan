classdef LHSintegratorAnisotropicStiffness < LHSintegrator
    
    properties (Access = private)
        geometry
        CAnisotropic
        Celas
        alphaDeg
        dofsInElem
    end

    methods (Access = public)
        
        function obj = LHSintegratorAnisotropicStiffness(cParams)
            obj.init(cParams);
            obj.initAnisotropicTensor(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end

        function LHS = compute(obj)
            obj.assemblyCMatrix();
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
        
    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dvolu = obj.mesh.computeDvolume(obj.quadrature);
            ngaus = obj.quadrature.ngaus;
            nelem = obj.mesh.nelem;
            ndpe  = obj.dim.ndofsElem;
            lhs = zeros(ndpe,ndpe,nelem);
            C   = obj.Celas;
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.computeBmat(igaus);
                for iel = 1:nelem
                    BmatEl = Bmat(:,:,iel);
                    dNdN = BmatEl'*C(:,:,iel)*BmatEl;
                    dV = dvolu(igaus,iel);
                    lhs(:,:,iel) = lhs(:,:,iel) + dNdN*dV;
                end
            end
        end
        
   end
    
   methods (Access = private)
       
       function initAnisotropicTensor(obj,cParams)
           CLocal = cParams.CAnisotropic;
           obj.alphaDeg = cParams.aniAlphaDeg;
           obj.CAnisotropic = obj.rotateAnisotropicMatrix(CLocal);
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

        function Bcomp = createBComputer(obj)
            s.dim          = obj.dim;
            s.geometry     = obj.geometry;
            s.globalConnec = obj.globalConnec;
            s.dofsInElem   = obj.dofsInElem;
            Bcomp = BMatrixComputer(s);
        end

        function CGlobal = rotateAnisotropicMatrix(obj,CLocal)
            R = [cosd(obj.alphaDeg),-sind(obj.alphaDeg)
                sind(obj.alphaDeg), cosd(obj.alphaDeg)];
            CGlobal = R*CLocal*R';
        end

        function assemblyCMatrix(obj)
            nelem = size(obj.mesh.connec,1);
            C = zeros(size(obj.CAnisotropic,1),size(obj.CAnisotropic,2),nelem);
            for i = 1:nelem
                C(:,:,i) = obj.CAnisotropic;
            end
            obj.Celas = C;
        end
   end

end