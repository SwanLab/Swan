classdef SymmetricPointComputer < handle
    
    properties (Access = private)
        tol
    end
    
    properties (Access = private)
        vector
        line
    end
    
    methods (Access = public)
        
        function obj = SymmetricPointComputer(cParams)
            obj.init(cParams)
        end
        
        function itIs = isNodeInLine(obj)
            dist = obj.computeDistanceToLine();
            itIs = dist < obj.tol;
        end
        
        function s = computeSymmetricVector(obj)
            pr = obj.computePRvector();
            p  = obj.vector;
            s = 2*pr + p;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.vector = cParams.vector;
            obj.line   = cParams.line;
            obj.tol    = 1e-13;
        end
        
        function pr = computePRvector(obj)
            pq = obj.computePQvector();
            rq = obj.computeRQvector(pq);
            qr = -rq;
            pr = pq + qr;
        end
        
        function dist = computeDistanceToLine(obj)
            pr = obj.computePRvector();
            dist = sqrt(pr(:,1).^2 + pr(:,2).^2);
        end
        
        function pq = computePQvector(obj)
            p = obj.vector;
            q = obj.line.point;
            ndim = size(p,2);
            pq = zeros(size(p));
            for idim = 1:ndim
                pq(:,idim) = q(idim) - p(:,idim);
            end
        end
        
        function rq = computeRQvector(obj,pq)
            u  = obj.line.vector;
            nnode = size(pq,1);
            uVect = repmat(u',nnode,1);
            proj = dot(pq,uVect,2);
            rq   = proj.*uVect;
        end
        
    end
    
end