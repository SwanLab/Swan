classdef MeshInterpolator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       fineMesh 
       coarseMesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MeshInterpolator()
            close all
            obj.createCoarseMesh()
            obj.createFunctions();
            obj.remeshDiscontinousField();
        end
        
    end
    
    methods (Access = private)
        
        
        function createCoarseMesh(obj)
            s.xmin = 0;
            s.xmax = 1;
            s.ymin = 0;
            s.ymax = 1;
            s.nInteriorPoints = 500;
            s.nBoundaryPoints = 15;
            uM = UnstructuredMeshCreator(s);
            m = uM.create();

% 
%             s.coord = [1 0; 0 1; 0 0; 1 1];
%             s.connec = [3 1 2; 1 4 2];
%             m = Mesh(s);

            obj.coarseMesh = m;
        end

        function mFine = remesh(obj,mC)
            %mC = obj.coarseMesh;
            mC.computeEdges();
            e = mC.edges;

            eMesh = mC.computeEdgeMesh();

            newCoord = eMesh.computeBaricenter()';
            allCoord = [mC.coord;newCoord];

            
            nnodes            = eMesh.nnodes;
            newNodes(:,1)     = (nnodes+1):(nnodes+eMesh.nelem);
            edgesWithNewNodes = 1:eMesh.nelem;
            edgesAsBefore     = setdiff(1:eMesh.nelem,edgesWithNewNodes);
            firstElems    = [eMesh.connec(edgesWithNewNodes,1),newNodes];
            secondElems   = [newNodes,eMesh.connec(edgesWithNewNodes,2)];
            newAllElems   = [firstElems,secondElems];
            newConecs     = reshape(newAllElems(:)',[],2);

            allConnecEdges = [eMesh.connec(edgesAsBefore,:),newConecs];                        

            sE.coord  = allCoord;
            sE.connec = allConnecEdges;
            sE.kFace  = eMesh.kFace;
            eMeshFine = Mesh(sE);
            figure()
            eMeshFine.plot()



            rCells = 1:mC.nelem;
            fCells = setdiff(1:mC.nelem,rCells);

            nV1 = mC.connec(rCells,1);
            nV2 = mC.connec(rCells,2);
            nV3 = mC.connec(rCells,3);
            
            nE1 = newNodes(e.edgesInElem(rCells,1),1);
            nE2 = newNodes(e.edgesInElem(rCells,2),1);
            nE3 = newNodes(e.edgesInElem(rCells,3),1);

            nConnec(:,1,:) = [nV1 nE1 nE3];
            nConnec(:,2,:) = [nE1 nV2 nE2];
            nConnec(:,3,:) = [nE1 nE2 nE3];
            nConnec(:,4,:) = [nE3 nE2 nV3];

            newConnec = reshape(nConnec,[],3);
            
            allConnec = [mC.connec(fCells,:),newConnec]; 

            figure()
            sF.connec = allConnec;
            sF.coord  = [mC.coord;newCoord];
            mFine = Mesh(sF);
            mFine.plot()
        end
        
        function createFunctions(obj,m)
            f   = @(x) x(:,1).*x(:,2);
            m   = obj.coarseMesh;
            fP1 = obj.createFeFunction(m,f);


            mFine     = obj.remesh(m);
            s.connec  = mFine.connec;
            s.type    = mFine.type;
            s.fValues = obj.interpolateFeFunction(f,m);
            cF = P1Function(s); 
            cF.plot(mFine) 

        end

        function fP1 = createFeFunction(obj,m,f)
            s.connec  = m.connec;
            s.type    = m.type;
            s.fValues = f(m.coord);
            fP1 = P1Function(s); 
            fP1.plot(m)
        end        

        function createFineFeFunction(obj)
            func = @(x) x(:,1).*x(:,2);            
            m = obj.coarseMesh;




            mD = obj.coarseMesh.createDiscontinuousMesh();
            mDFine = obj.remesh(mD);
        end




        function allF = interpolateFeFunction(obj,f,m)
            me = m.computeEdgeMesh();  
            fe = obj.createFeFunction(me,f); 
            q = Quadrature.set(m.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fNewValues = squeeze(fe.evaluate(xV));                
            allF = [fe.fValues;fNewValues];
        end

        function remeshDiscontinousField(obj)
            mC = obj.coarseMesh;
            m = mC.createDiscontinuousMesh();
            figure()
            m.plot()

            f = obj.createFunction(mC);

            mFine = obj.fineMesh;
            mFine.createDiscontinuousMesh();
        
            e = mC.edges;
            m.computeEdges;
            eD = m.edges;
            connec = m.connec;
            for iEdge = 1:3
                nodes = e.localNodeByEdgeByElem(1,iEdge,:);
                nodeA = nodes(:,1);
                nodeB = nodes(:,2);


            end

        end

        function fD = createFunction(obj,m)
            %m = m.createDiscontinuousMesh();
            coord = m.computeBaricenter()';
            ls = coord(:,2)-(1-coord(:,1));
            chi = zeros(size(ls));
            chi(ls<0) = 1;            
            %= heaviside(ls);

            

            s.fValues = chi;%reshape(chi,[],3); 
            s.connec  = m.connec; 
            s.type    = m.type;
            f = P0Function(s);

            s.mesh = m;
            s.connec = m.connec;
            p = Projector_toP1Discontinuous(s);
            fD = p.project(f);
            fD.plot(m)
            end


        
    end
    
end