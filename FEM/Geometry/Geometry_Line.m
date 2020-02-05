classdef Geometry_Line < Geometry
    
    properties (GetAccess = public, SetAccess = private)
        xGauss
        dvolu
    end
    
    properties (Access = private)
        interpolationGeometry
        interpolationVariable  
        quadrature
        nElem
        coordElem        
    end
    
    methods (Access = public)
        
        function obj = Geometry_Line(cParams)
            obj.init(cParams)
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.computeElementCoordinates();
            obj.computeGaussPointsPosition();
            obj.computeDvolu(igaus);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.mesh.nelem;
            obj.interpolationGeometry = Interpolation.create(cParams.mesh,'LINEAR');
        end
        
        function initGeometry(obj,interpV,quad)
            obj.interpolationVariable = interpV;
            obj.quadrature = quad;
            obj.computeShapeFunctions();
        end
        
        function computeShapeFunctions(obj)
            xpg = obj.quadrature.posgp;
            obj.interpolationVariable.computeShapeDeriv(xpg)
            obj.interpolationGeometry.computeShapeDeriv(xpg);
        end
        
        function computeElementCoordinates(obj)
            nNode  = obj.interpolationGeometry.nnode;
            nDime  = obj.interpolationGeometry.ndime;
            coord  = obj.interpolationGeometry.xpoints;
            connec = obj.interpolationGeometry.T;
            coordE = zeros(nNode,nDime,obj.nElem);
            coord  = coord';
            for inode = 1:nNode
                nodes = connec(:,inode);
                coordNodes = coord(:,nodes);
                coordE(inode,:,:) = coordNodes;
            end
            obj.coordElem = coordE;
        end 
        
        function computeGaussPointsPosition(obj)
            nNode  = obj.interpolationGeometry.nnode;
            nDime  = obj.interpolationGeometry.ndime;
            shapes = obj.interpolationGeometry.shape;
            nGaus  = obj.quadrature.ngaus;
            xGaus = zeros(nGaus,nDime,obj.nElem);
            for kNode = 1:nNode
                shapeKJ(:,1) = shapes(kNode,:)';
                xKJ = obj.coordElem(kNode,:,:);
                xG = bsxfun(@times,shapeKJ,xKJ);
                xGaus = xGaus + xG;
            end
            obj.xGauss = permute(xGaus,[2 1 3]);
        end     
        
        function computeDvolu(obj)
            nGaus  = obj.quadrature.ngaus;
            nDime  = obj.interpolationGeometry.ndime;            
            drDtxi = zeros(nGaus,obj.nElem);
            xp     = permute(obj.coordElem,[1 3 2]);
            deriv  = obj.interpolationGeometry.deriv(1,:,:);
            dShapes = permute(deriv,[3 2 1]);
            for idime = 1:nDime
                x(:,1) = xp(:,:,idime);
                dxDtxi = dShapes*x;
                drDtxi = drDtxi + (dxDtxi).^2;
            end
            w(:,1) = obj.quadrature.weigp;
            obj.dvolu = bsxfun(@times,w,sqrt(drDtxi));
        end        
        
    end
    
    
    
end