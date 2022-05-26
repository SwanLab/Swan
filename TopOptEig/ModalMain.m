classdef ModalMain < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        dim
        material
        interpolation
        coord 
        connec
        mesh
        freeNodes
        stiffnessMatrixComputer
        massMatrixComputer
    end
    
    methods (Access = public)
        
        function obj = ModalMain()
            obj.init()
            obj.createMesh();
            dE = obj.createDensity();
            mI = obj.createMaterialInterpolation(dE);
            obj.createMaterial(mI);
            obj.createLHS();
            obj.solveEigModes();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            x = 1:10;
            y = 1:50;
            [X,Y] = meshgrid(x,y);
            i = 1;
             for iRow = 1:size(X,1)
                for iColum = 1: size(X,2)
                xc = X(iRow,iColum);
                yc = Y(iRow,iColum);
                obj.coord(i,:) = [xc,yc];
                i = i+1;
                end
             end
            DT = delaunay(X,Y);
            obj.connec = DT;
            triplot(DT,X,Y)
        end
        
        function createMesh(obj)
            s.coord  = obj.coord();
            s.connec = obj.connec();
            s.type = 'LINE';
            m = Mesh(s);
            obj.mesh = m;
        end

        function mI = createMaterialInterpolation(obj,dE)
            sD.ngaus = 1;
            sD.mesh  = obj.mesh;
            sD.type  = 'Vector';
            sD.ndimf = 2;
            sD.fieldName = 'Disp';
            d = DimensionVariables.create(sD);
            d.compute();
            obj.dim = d;
            s.dim = '2D';
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.constitutiveProperties.rho_plus = 1;
            s.constitutiveProperties.rho_minus = 1e-3;
            s.constitutiveProperties.E_plus = 1;
            s.constitutiveProperties.E_minus = 1e-3;
            s.constitutiveProperties.nu_minus = 1/3;
            s.constitutiveProperties.nu_plus = 1/3;
            s.nElem = obj.mesh.nelem;
            m = MaterialInterpolation.create(s);
            mI = m.computeMatProp(dE);
            obj.interpolation = mI;
        end

        function dE = createDensity(obj)
            m = obj.mesh;
            s.type = 'Density';
            s.mesh = m;
            s.value = [];
            s.initialCase = 'circleInclusion';
            s.creatorSettings.type = 'FromLevelSet';
            s.creatorSettings.fracRadius = 0.5;
            d = Density(s);
            
            s.mesh = m;
            s.field = d.value;
            p = NodalFieldPlotter(s);
            p.plot()
            
            sF.connec = m.connec;
            sF.type   = m.type;
            sF.fNodes = d.value;
            d = FeFunction(sF);
            dE = d.computeValueInCenterElement();
        end

        function createMaterial(obj,mI)
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mI.kappa;
            s.mu    = mI.mu;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end

        function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.ndof-1,d.ndof]);
            nodes = 1:d.ndof;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
        end  

        function createLHS(obj)
            sS.type = 'ElasticStiffnessMatrix';
            sS.dim = obj.dim;
            sS.mesh = obj.mesh;
            sS.globalConnec   = obj.mesh.connec;
            sS.freeNodes      = obj.freeNodes;
            sS.material       = obj.material;
            sS.interpolation  = obj.interpolation;
            obj.stiffnessMatrixComputer = LHSintegrator.create(sS);

%             sM.type         = 'MassMatrix';
%             sM.dim          = obj.dim;
%             sM.mesh         = obj.mesh;
%             sM.globalConnec = obj.mesh.connec;
%             sM.freeNodes    = obj.freeNodes;
%             sM.quadType     = 'LINE';
%             obj.massMatrixComputer = LHSintegrator.create(sM);
        end

        function solveEigModes(obj)
%             s.mesh = obj.mesh;
%             s.dim  = obj.dim;
%             s.stiffnessMatrix = obj.stiffnessMatrix;
%             s.massMatrix = obj.massMatrix;
            K = obj.stiffnessMatrixComputer.compute();
           % M = obj.massMatrixComputer.compute();
            Kfree = obj.provideFreeMatrix(K);
            Mfree = obj.provideFreeMatrix(M);
            [v,d] = eigs(Mfree,Kfree,2,'SM');
            obj.V  = v;
            obj.D  = d;
        end

        function MatrixFree = provideFreeMatrix(obj,Matrix)
            free = obj.free;
            MatrixFree = Matrix(free,free);
        end
        
    end
    
end