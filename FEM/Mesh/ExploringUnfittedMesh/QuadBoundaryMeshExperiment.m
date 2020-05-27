classdef QuadBoundaryMeshExperiment < MeshAndGaussOfUnfittedQuadMesh
    
    methods (Access = public)
        
        function obj = QuadBoundaryMeshExperiment()
            
            obj.init();
            close all
            
%             coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
%             connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
            ls = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
            
            nameCase = 'cutMeshProvisionalQuadBoundary';
            s.coord = obj.backgroundMesh.coord;
            s.connec = obj.backgroundMesh.connec;
            m = Mesh().create(s);
            
            figure()
            s.meshBackground = m;
            s.unfittedType = 'INTERIOR';
            uMesh = UnfittedMesh(s);
            uMesh.compute(ls);
            uMesh.plot();
            
            [connecFull,connecCut,cutElems] = obj.computeConnecCutAndFull(obj.levelSet,obj.backgroundMesh.connec);
            
            sM.coord = obj.backgroundMesh.coord;
            sM.connec = connecCut;
            backgroundCutMesh = Mesh().create(sM);
            
            s.backgroundMesh = backgroundCutMesh;
            s.cutElems = cutElems;
            s.levelSet  = ls;
            s.lastNode = max(obj.backgroundMesh.connec(:));
            cutMesh = CutMeshProvisionalQuadrilater(s);
            cutMesh.compute();
            
            
            quad = Quadrature.set(uMesh.innerCutMesh.geometryType);
            quad.computeQuadrature('QUADRATIC2');
            m = uMesh.innerCutMesh.cutMeshOfSubCellLocal;
            xG2 = m.computeXgauss(quad.posgp);
            %xG3 = reshape(xG2,2,[])';
            
            s.coord = m.coord;
            conC = obj.backgroundMesh.connec(uMesh.innerCutMesh.cellContainingSubcell,:);
            s.connec = conC;
            m = Mesh().create(s);
            
            xC = uMesh.innerCutMesh.cutMeshOfSubCellGlobal.computeXgauss(xG2);
            %figure(10)
            hold on
            plot(xC(1,:),xC(2,:),'*','MarkerSize',10)
            
            
            bMesh = cutMesh.computeBoundaryMesh();
            
            quad = Quadrature.set(bMesh.geometryType);
            quad.computeQuadrature('QUADRATIC');
            
            
            xCutB = bMesh.computeXgauss(quad.posgp);
            plot(xCutB(1,:),xCutB(2,:),'g*','MarkerSize',10)
            
            patch('vertices',bMesh.coord,'faces',bMesh.connec,...
                'edgecolor',[0 0 1], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat','LineWidth',2)
            axis('equal');
            
            d = load(nameCase);
            %
            error(1) = norm(d.boundaryMesh.connec(:) - bMesh.connec(:));
            error(2) = norm(d.boundaryMesh.coord(:)  - bMesh.coord(:));
            error(3) = norm(d.interiorMesh.connec(:)    - uMesh.innerMesh.connec(:));
            error(4) = norm(d.interiorMesh.coord(:)     - uMesh.innerMesh.coord(:));
            error(5) = norm(d.interiorCutMesh.connec(:) - uMesh.innerCutMesh.connec(:));
            error(6) = norm(d.interiorCutMesh.coord(:)  - uMesh.innerCutMesh.coord(:));
            error(7) = norm(d.xGaussInteriorCutMesh(:) - xC(:));
            error(8) = norm(d.xGaussBoundaryMesh(:) - xCutB(:));
            
            errorNorm = norm(error)
            
            
            
        end
        

                
        
    end
    
end
