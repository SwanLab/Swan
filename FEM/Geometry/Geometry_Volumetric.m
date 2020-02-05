classdef Geometry_Volumetric < Geometry
    
    properties (GetAccess = public, SetAccess = private)
        cartd
        dvolu
    end
    
    properties (Access = private)
        mesh
        interpolationVariable
        quadrature                
        matrixInverter
        jacobian
        detJ
    end
    
    methods (Access = public)
        
        function obj = Geometry_Volumetric(cParams)
            obj.init(cParams)
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            for igaus = 1:obj.quadrature.ngaus
                obj.computeJacobian(igaus);
                obj.computeJacobianDeterminant(igaus);
                obj.computeDvolu(igaus);
                obj.computeCartesianDerivatives(igaus);
            end
        end
        
    end
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function initGeometry(obj,interpV,quad)
            obj.interpolationVariable = interpV;
            obj.quadrature = quad;
            obj.computeShapeFunctions();
            obj.initDvolu();
            obj.initDetJ();
            obj.initCartD();
            obj.matrixInverter = MatrixVectorizedInverter();
        end
        
        function computeShapeFunctions(obj)
            xpg = obj.quadrature.posgp;
            obj.interpolationVariable.computeShapeDeriv(xpg)
            obj.mesh.interpolation.computeShapeDeriv(xpg)            
        end
        
        function initDvolu(obj)
            nGaus     = obj.quadrature.ngaus;
            obj.dvolu = zeros(obj.mesh.nelem,nGaus);
        end
        
        function initDetJ(obj)
            nGaus    = obj.quadrature.ngaus;
            obj.detJ = zeros(obj.mesh.nelem,nGaus);
        end
        
        function initCartD(obj)
            nDime = obj.interpolationVariable.ndime;
            nNode = obj.interpolationVariable.nnode;
            nGaus = obj.quadrature.ngaus;
            obj.cartd  = zeros(nDime,nNode,obj.mesh.nelem,nGaus);
        end
          
        function computeJacobian(obj,igaus)
            nDime   = obj.mesh.ndim;
            nNode   = obj.mesh.nnode;
            nElem   = obj.mesh.nelem;
            dShapes = obj.mesh.interpolation.deriv(:,:,igaus);
            jac = zeros(nDime,nDime,nElem);
            for kNode = 1:nNode
                dShapeIK = dShapes(:,kNode);
                xKJ      = obj.mesh.coordElem(kNode,:,:);
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
            nElem   = obj.mesh.nelem;
            nNode   = obj.interpolationVariable.nnode;
            nDime   = obj.interpolationVariable.ndime;
            dShapes = obj.interpolationVariable.deriv(:,:,igaus);
            invJ     = obj.computeInvJacobian();
            dShapeDx = zeros(nDime,nNode,nElem);
            for jDime = 1:nDime
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