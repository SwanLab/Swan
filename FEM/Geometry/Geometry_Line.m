classdef Geometry_Line < Geometry
    
    properties (Access = public)
        normalVector        
    end
    
    properties (Access = private)
        drDtxi
        tangentVector
    end
    
    methods (Access = public)
        
        function obj = Geometry_Line(cParams)
            obj.permutation = [2 3 1];
            obj.init(cParams);
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.computeDvolu();
        end
        
    end
    
    methods (Access = private)
        
        function computeDrDtxi(obj)
            nGaus  = obj.quadrature.ngaus;
            nDime  = obj.mesh.ndim;
            drdtxi = zeros(nGaus,obj.mesh.nelem,nDime);
            xp = obj.coordElem;
            deriv  = obj.mesh.interpolation.deriv(1,:,:);
            dShapes = permute(deriv,[3 2 1]);
            for idime = 1:nDime
                x      = xp(:,:,idime);
                drdtxi(:,:,idime) = dShapes*x;
            end
            obj.drDtxi = drdtxi;
        end
        
        function computeTangentVector(obj)
            obj.computeDrDtxi();
            drdtxi = obj.drDtxi;
            drdtxiNorm = obj.computeVectorNorm(drdtxi);
            t = drdtxi./drdtxiNorm;
            obj.tangentVector = t;
        end
        
        function computeNormalVector(obj)
            nDime  = obj.mesh.ndim;
            switch nDime
                case 2
                    error('Not done, Call second derivative');
                case 3
                    obj.computeTangentVector();
                    t = obj.tangentVector;
                    n = zeros(size(t));
                    n(:,:,1) = -t(:,:,2);
                    n(:,:,2) =  t(:,:,1);
                    obj.normalVector = n;
            end
        end
        
        function computeDvolu(obj)
            obj.computeDrDtxi();
            drdtxiNorm = obj.computeVectorNorm(obj.drDtxi);
            w(:,1) = obj.quadrature.weigp;
            dv =  bsxfun(@times,w,drdtxiNorm);
            obj.dvolu = dv';
        end
        
    end
    
    methods (Access = private, Static)
        
        function nv = computeVectorNorm(v)
            nDime  = size(v,3);
            nv = zeros(size(v,1),size(v,2));
            for idime = 1:nDime
                ti    = v(:,:,idime);
                nv = nv + (ti).^2;
            end
            nv = sqrt(nv);
        end
        
    end
    
end