classdef AbstractMesh < handle
    
%     properties (Access = public)
%         unfittedType        
%     end
    
    
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        
        nelem
        ndim
        
        geometryType
        
        coordElem
        interpolation
    end
    
    methods (Access = public)
       
       function xGauss = computeXgauss(obj,quad)
            xpg = quad.posgp;
            obj.interpolation.computeShapeDeriv(xpg);           
            nNode  = obj.nnode;
            nDime  = obj.ndim;
            shapes = obj.interpolation.shape;
            nGaus  = quad.ngaus;
            xGauss = zeros(nGaus,nDime,obj.nelem);
            for kNode = 1:nNode
                shapeKJ(:,1) = shapes(kNode,:)';
                xKJ = obj.coordElem(kNode,:,:);
                xG = bsxfun(@times,shapeKJ,xKJ);
                xGauss = xGauss + xG;
            end
            xGauss = permute(xGauss,[2 1 3]);
       end        
       
       function dvolume = computeDvolume(obj,quad)
            s.mesh = obj;
            g = Geometry.create(s);
            g.computeGeometry(quad,obj.interpolation);
            dvolume = g.dvolu;
            dvolume = dvolume';
        end       
                
    end
    
    methods (Access = protected)
        
        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj,'LINEAR');            
        end        
        
        function computeElementCoordinates(obj)
            nNode  = obj.nnode;
            nDime  = obj.ndim;
            coordE = zeros(nNode,nDime,obj.nelem);
            coords = obj.coord';
            for inode = 1:nNode
                nodes = obj.connec(:,inode);
                coordNodes = coords(:,nodes);
                coordE(inode,:,:) = coordNodes;
            end
            obj.coordElem = coordE;
        end        
        
        
    end
    
end