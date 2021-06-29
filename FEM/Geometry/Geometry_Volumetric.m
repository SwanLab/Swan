classdef Geometry_Volumetric < Geometry
    
    properties (GetAccess = public, SetAccess = private)
        cartd
    end
    
    properties (Access = private)              
        matrixInverter
        jacobian
        detJ
    end
    
    methods (Access = public)
        
        function obj = Geometry_Volumetric(cParams)
            obj.permutation = [2 1 3];                        
            obj.init(cParams);
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.initVariables();
            obj.matrixInverter = MatrixVectorizedInverter();            
            for igaus = 1:obj.quadrature.ngaus
                obj.computeJacobian(igaus);
                obj.computeJacobianDeterminant(igaus);
                obj.computeDvolu(igaus);
                obj.computeCartesianDerivatives(igaus);
            end
        end
        
    end
    
    methods (Access = private)
    
        function initVariables(obj)
            nDime = obj.interpolationVariable.ndime;
            nNode = obj.interpolationVariable.nnode;
            nGaus = obj.quadrature.ngaus;
            obj.cartd = zeros(nDime,nNode,obj.mesh.nelem,nGaus);
            obj.detJ  = zeros(obj.mesh.nelem,nGaus);
            obj.dvolu = zeros(obj.mesh.nelem,nGaus);  
        end
          
        function computeJacobian(obj,igaus)
            nDime   = obj.mesh.ndim;
            nNode   = obj.mesh.nnode;
            nElem   = obj.mesh.nelem;
            dShapes = obj.mesh.interpolation.deriv(:,:,igaus);
            jac = zeros(nDime,nDime,nElem);           
            for kNode = 1:nNode
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