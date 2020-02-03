classdef Geometry < handle
    
    properties (GetAccess = public, SetAccess = public)
        cartd
        cart_pos_gp
        dvolu        
        type
        nfields

    end
    
    properties (Access = private)
        matrixInverter
        interpolation                                     
        interpolationVariable
        quadrature
        ngaus
        nnode
        ndime
        nelem
        coordElem
        jacobian
        detJ
    end
    
    methods (Access = public)
        
        function obj = Geometry(mesh,order)
            obj.type = mesh.geometryType;   
            obj.nfields = 1;            
            %obj.interpolation = Interpolation.create(mesh,order); %!!!!!!!!!!!
            obj.interpolation = Interpolation.create(mesh,'LINEAR');            
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.init(interpV,quad);
            obj.computeElementCoordinates();
            obj.computeGaussPointsPosition();
            for igaus = 1:obj.ngaus
                obj.computeJacobian(igaus);
                obj.computeJacobianDeterminant(igaus);
                obj.computeDvolu(igaus);
                obj.computeCartesianDerivatives(igaus);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,interpV,quad)
            obj.interpolationVariable = interpV;
            obj.quadrature = quad;           
            obj.ndime  = obj.interpolationVariable.ndime;
            obj.nnode  = obj.interpolationVariable.nnode;
            obj.ngaus  = obj.quadrature.ngaus;
            obj.nelem  = obj.interpolationVariable.nelem;
            obj.dvolu  = zeros(obj.nelem,obj.ngaus);
            obj.detJ   = zeros(obj.nelem,obj.ngaus);
            obj.cartd  = zeros(obj.ndime,obj.nnode,obj.nelem,obj.ngaus);
            obj.matrixInverter = MatrixVectorizedInverter();            
        end
        
        function computeGaussPointsPosition(obj)
            nNode = obj.interpolationVariable.nnode;
            shapes = obj.interpolationVariable.shape;
            xGauss = zeros(obj.ngaus,obj.ndime,obj.nelem);
            for kNode = 1:nNode 
                shapeKJ(:,1) = shapes(kNode,:)';
                xKJ = obj.coordElem(kNode,:,:);
                xG = bsxfun(@times,shapeKJ,xKJ);
                xGauss = xGauss + xG;                                                
            end
            obj.cart_pos_gp = permute(xGauss,[2 1 3]);
        end
        
        function computeElementCoordinates(obj)
            coord  = obj.interpolationVariable.xpoints;
            connec = obj.interpolationVariable.T;
            coordE = zeros(obj.nnode,obj.ndime,obj.nelem);
            nDime  = size(coord,2);
            for idime = 1:nDime
                coordDim = coord(:,idime);
                for inode = 1:obj.nnode
                    coordE(inode,idime,:) = coordDim(connec(:,inode));
                end
            end
            obj.coordElem = coordE;
        end
        
        function computeJacobian(obj,igaus)
            jac = zeros(obj.ndime,obj.ndime,obj.nelem);
            dShapes = obj.interpolationVariable.deriv(:,:,igaus);
            for kNode = 1:obj.nnode
                dShapeIK = dShapes(:,kNode);
                xKJ      = obj.coordElem(kNode,:,:);
                jacIJ    = bsxfun(@times, dShapeIK, xKJ);
                jac = jac + jacIJ;
            end
            obj.jacobian = jac;
        end
        
        function computeJacobianDeterminant(obj,igaus)
            J = obj.jacobian;
            obj.detJ(:,igaus) = obj.matrixInverter.computeDeterminant(J);
        end
        
        function computeDvolu(obj,igaus)
            w = obj.quadrature.weigp;
            obj.dvolu(:,igaus) = w(igaus)*obj.detJ(:,igaus);
        end
        
        function computeCartesianDerivatives(obj,igaus)
            invJ     = obj.computeInvJacobian();
            dShapes  = obj.interpolationVariable.deriv(:,:,igaus);
            dShapeDx = zeros(obj.ndime,obj.nnode,obj.nelem);
            for jDime = 1:obj.ndime
                invJ_JI   = invJ(:,jDime,:);
                dShape_KJ = dShapes(jDime,:);
                dSDx_KI   = bsxfun(@times, invJ_JI,dShape_KJ);
                dShapeDx  = dShapeDx + dSDx_KI;
            end
            obj.cartd(:,:,:,igaus) = dShapeDx;
        end
        
        function invJ = computeInvJacobian(obj)
            jac = obj.jacobian;
            invJ = obj.matrixInverter.computeInverse(jac);
        end
        
    end
    
end