classdef SurfaceMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Surface';
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = SurfaceMesh(cParams)
            obj = obj@Mesh(cParams);            
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            switch obj.ndim
                case 3 % 2D Surface in 3D space
                    J = obj.computeJacobian(xV);
                    n = obj.computeNormalVectors(J);
                    detJ = squeeze(pagenorm(n));
                %case 2 % 2D Surface in 2D space
                %    J = obj.computeJacobian(xV);
                %    detJ = MatrixVectorizedInverter.computeDeterminant(J);
            end
        end

        function n = computeNormals(obj, xV)
            J = obj.computeJacobian(xV);
            n = obj.computeNormalVectors(J);
        end

        function n = getNormals(obj) 
            quad = Quadrature.create(obj,0);
            n = obj.computeNormals(quad.posgp);
        end

        function plotSolidColor(obj,color)
            faceColor = color;
            faceAlpha = 1;
            edgeAlpha = 0;
            obj.plotSpecific(faceColor,faceAlpha,edgeAlpha)
        end    

        function plot(obj) 
            faceColor = "red";
            faceAlpha = 0.3;
            edgeAlpha = 0.5;
            obj.plotSpecific(faceColor,faceAlpha,edgeAlpha)
        end

        function m = provideExtrudedMesh(obj, height)
            s.unfittedMesh = obj;
            s.height       = height;
            me = MeshExtruder(s);
            m = me.extrude();
        end
        
    end
    
    methods (Access = private)

        function plotSpecific(obj,faceColor,faceAlpha,edgeAlpha)
            if size(obj.connec,2) == 3 && size(obj.coord,2) == 3
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                z = obj.coord(:,3);
                p = trisurf(obj.connec,x,y,z);
                p.FaceColor = [1 0 0];
                p.FaceAlpha = 1;
                p.EdgeColor = 'none';
                hold on
            else
                p = patch('vertices',obj.coord,'faces',obj.connec);
                p.EdgeAlpha = edgeAlpha;
                p.EdgeLighting = 'flat';
                p.FaceColor = faceColor;%[167,238,237]/265; 'green';'red';%
                p.FaceLighting = 'flat';
                p.FaceAlpha = faceAlpha;
                p.LineWidth = 1.5;
                axis('equal');
                hold on
            end
        end        
        
        function normalVector = computeNormalVectors(obj,J)
            nDimGlo = size(J,2);
            nPoints = size(J,3);
            nElem = size(J,4);

            normalVector = zeros(1,nDimGlo,nPoints,nElem);
            % DxDxi  = squeeze(J(1,1,:,:))';
            % DxDeta = squeeze(J(2,1,:,:))';
            % DyDxi  = squeeze(J(1,2,:,:))';
            % DyDeta = squeeze(J(2,2,:,:))';
            % DzDxi  = squeeze(J(1,3,:,:))';
            % DzDeta = squeeze(J(2,3,:,:))';

            DxDxi  = squeezeParticular(J(1,1,:,:), [1 2]);
            DxDeta = squeezeParticular(J(2,1,:,:), [1 2]);
            DyDxi  = squeezeParticular(J(1,2,:,:), [1 2]);
            DyDeta = squeezeParticular(J(2,2,:,:), [1 2]);
            DzDxi  = squeezeParticular(J(1,3,:,:), [1 2]);
            DzDeta = squeezeParticular(J(2,3,:,:), [1 2]);
            normalVector(:,1,:,:) = DyDxi.*DzDeta - DzDxi.*DyDeta;
            normalVector(:,2,:,:) = DzDxi.*DxDeta - DxDxi.*DzDeta;
            normalVector(:,3,:,:) = DxDxi.*DyDeta - DyDxi.*DxDeta;
        end
        
    end
    
end