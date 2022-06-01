classdef ModalMain < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        dim
        material
        materialInterpolation
        interpolation
        coord 
        connec
        mesh
        freeDOFs
        stiffnessMatrixComputer
        massMatrixComputer
        V
        D
        mode1
        mode2
    end
    
    methods (Access = public)
        
        function obj = ModalMain()
            obj.init()
            obj.createMesh();
            dE = obj.createDensity();
            mI = obj.createMaterialInterpolation(dE);
            obj.createMaterial(mI);
            obj.createBoundaryConditions();
            obj.createLHS();
            obj.solveEigModes();
            obj.plotEigModes();
            obj.printInGiD();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            x1        = linspace(0,3,20);
            x2        = linspace(0,2,20);
            [X,Y]     = meshgrid(x1,x2);
            Z         = zeros(size(X));
            xmesh = X;
            ymesh = Y;
            zmesh = Z;

            [F,V]    = mesh2tri(xmesh,ymesh,zmesh,'f');
            obj.coord  = V(:,1:2);
            obj.connec = F;
        end
        
        function createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
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
            obj.materialInterpolation = mI;
        end

        function dE = createDensity(obj)
            m = obj.mesh;
            s.type = 'Density';
            s.mesh = m;
            s.value = [];
            s.initialCase = 'full';
            s.creatorSettings.type = 'FromLevelSet';
            s.creatorSettings.fracRadius = 0.5;
            d = Density(s);
            
%             s.mesh = m;
%             s.field = d.value;
%             p = NodalFieldPlotter(s);
%             p.plot()
            
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
            FixNod = obj.computeFixedNodes();
            FixDof = obj.computeFixedDOFs(FixNod); 
            dofs = 1:d.ndofs;
            free  = setdiff(dofs,FixDof);
            obj.freeDOFs = free;
        end  

        function FN = computeFixedNodes(obj)
            coorX = obj.coord(:,1);
            tol = 0.001;
            minX = min(coorX) + tol;
            maxX = max(coorX) - tol;
            FN1 = find(coorX < minX); 
            FN2 = find(coorX > maxX);
            FN = [FN1 ;FN2];
        end

        function FixDof = computeFixedDOFs(obj, FixNod)
            d = obj.dim;
            nDOFn = d.ndofs/d.nnodes;
            FD = zeros(nDOFn*length(FixNod),1);
            for i = 1: length(FixNod)
               FD(2*i-1,1) = FixNod(i)*2-1;
               FD(2*i,1) = FixNod(i)*2;
            end
            FixDof = FD';
        end

        function createLHS(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATIC');
            sM.quadrature = quad;
            sM.quadType     = 'QUADRATIC';
            int.computeShapeDeriv(quad.posgp);
            sM.interpolation  = int;
            sS.type = 'ElasticStiffnessMatrix';
            sS.dim = obj.dim;
            sS.mesh = obj.mesh;
            sS.globalConnec   = obj.mesh.connec;
            sS.freeNodes      = obj.freeDOFs;
            sS.material       = obj.material;
            sS.interpolation  = int;
            obj.stiffnessMatrixComputer = LHSintegrator.create(sS);

            int = Interpolation.create(obj.mesh,'LINEAR');
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATICMASS');
            sM.quadrature = quad;
            sM.quadType     = 'QUADRATICMASS';
            int.computeShapeDeriv(quad.posgp);
            sM.interpolation  = int;
            sM.type         = 'MassMatrix';
            sM.dim          = obj.dim;
            sM.mesh         = obj.mesh;
            sM.globalConnec = obj.mesh.connec;
            sM.freeNodes    = obj.freeDOFs;
            obj.massMatrixComputer = LHSintegrator.create(sM);
        end

        function solveEigModes(obj)
            K = obj.stiffnessMatrixComputer.compute();
            M = obj.massMatrixComputer.compute();
            Kfree = obj.provideFreeMatrix(K);
            Mfree = obj.provideFreeMatrix(M);
            [v,d] = eigs(Kfree,Mfree,2,'SM');
            obj.V  = v;
            obj.D  = d;
            obj.computeBucklingModes();
            obj.computeDisplacement(M);
        end

        function MatrixFree = provideFreeMatrix(obj,Matrix)
            free = obj.freeDOFs;
            MatrixFree = Matrix(free,free);
        end
        
        function [Mode1,Mode2] = computeBucklingModes(obj)
            nnode = obj.mesh.nnodes;
            Mode1=zeros(nnode,2);
            Mode2=zeros(nnode,2);
            v1 = obj.V(:,1);
            v2 = obj.V(:,2);
            for i=1 : nnode-40
            Mode1(i+40,:) = [v1(2*i-1) v1(2*i)];
            Mode2(i+40,:) = [v2(2*i-1) v2(2*i)];
            end
            obj.mode1 = Mode1; 
            obj.mode2 = Mode2; 
        end 

        function computeDisplacement(obj,M)
            coor = obj.mesh.coord;
            mod1entero = zeros(obj.dim.ndofs,1);
            mod1entero(obj.freeDOFs') = obj.V(:,1);
            t_m = mod1entero.';
            d = zeros(size(obj.dim.ndofs,1),1);
            for i = 1: size(coor,1)
                FD(2*i-1,1) = coor(i,1);
                FD(2*i,1) = coor(i,2);
            end
            disp = t_m*M*FD;
            disp = mod1entero*disp.*FD;
        end

        function plotEigModes(obj)
            mod1 = obj.mode1;
            coord = obj.mesh.coord;
            subplot(2,1,1); plot( mod1(:,1), mod1(:,2));
            grid on
            grid minor
            title('First Mode X','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,1,2); plot(coord(:,2), mod1(:,2));
            grid on
            grid minor
            title('First Mode Y')
        end

        function printInGiD(obj)  %% no puedo seguir
            fileName = 'Eig';
            m = obj.mesh;
            quad = Quadrature.set(m.type);
            quad.computeQuadrature('LINEAR');
            dI = obj.createPostProcessDataBase(m,fileName);
            dI.ndim   = 2;
            dI.pdim   = '2D';
            dI.ptype  = 'ELASTIC';
            dI.name = '';
            dI.printers = 'DensityGauss'; % 'HomogenizedTensor'
            p = Postprocess('VectorField',dI);
            %f.u = obj.mode1;
            dP.fields.u = obj.mode1;
            dP.quad   = quad;
            iter = 0;
            p.print(iter,dP);
        end

        function d = createPostProcessDataBase(obj,mesh,fileName)
            dI.mesh    = mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
        end



    end

end