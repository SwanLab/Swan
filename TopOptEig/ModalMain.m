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
            s.initialCase = 'circle';
            s.creatorSettings.type = 'FromLevelSet';
            s.creatorSettings.fracRadius = 0.5;
            d = Density(s);
            
            %s.mesh = m;
            %s.field = d.value;
            %p = NodalFieldPlotter(s);
            %p.plot()
            
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
            sS.interpolation  = obj.interpolation;
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
            subplot(2,2,2); plot( mod1(:,1), mod1(:,2));
            grid on
            grid minor
            title('First Mode X','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(coord(:,2), mod1(:,2));
            grid on
            grid minor
            title('First Mode Y')
        end

        %         function printInGiD(obj)  %% no puedo seguir
        %             mesh = obj.mesh;
        %             s.material = obj.material;
        %             s.type = 'ELASTIC';
        %             s.scale = 'MACRO';
        %             fem = FEM.create(s);
        %             fem.setC(mat.C);
        %             fem.computeChomog();
        %             quad = Quadrature.set(mesh.type);
        %             quad.computeQuadrature('LINEAR');
        %             dI = createPostProcessDataBase(mesh,'');
        %             dI.ndim   = 2;
        %             dI.pdim   = '2D';
        %             dI.ptype  = 'ELASTIC';
        %             dI.name = '';
        %             dI.printers = 'DensityGauss'; % 'HomogenizedTensor'
        %             p = Postprocess('NumericalHomogenizer',dI);
        %             dP.fields{1} = fem.variables2print;
        %             dP.fields{2} = Modes;
        %             %dP.fields{2} = density;
        %             dP.quad   = quad;
        %             iter = 0;
        %             p.print(iter,dP);
        %         end

        function printInGiD(obj)
            coor = obj.coord;
            conn = obj.connec;
            MODES = obj.V;
            posgp = 1;
            DOFl = obj.freeDOFs;
            TypeElement = 'TRIANGLE';
            NameFileMesh = 'PostModes';
            NameFile_msh = ['GIDPOST/','MODES','_',NameFileMesh,'.msh'] ;
            NameFile_res= ['GIDPOST/','MODES','_',NameFileMesh,'.res'] ;
            MaterialType = ones(size(conn,1),1) ; 
            obj.GidMesh2DFE(NameFile_msh,coor,conn,MaterialType,TypeElement);
            MODESplot = zeros(size(coor,1)*size(coor,2),size(MODES,2)) ;
            MODESplot(DOFl,:) = MODES;
            obj.GidResults2DFE_modes(NameFile_res,coor,TypeElement,MODESplot,posgp);
            cddd = cd ;
            NAMEFILEOPEN =  [cddd,'/',NameFile_res] ;
            disp('open GID FILE:')
            disp(NAMEFILEOPEN)
        end

        function GidMesh2DFE(obj,NameFile,COOR,CONNECT,MaterialType,TypeElement)
            NAMEPROJ = 'MODES';
            numer_nodes = 1:size(COOR,1);
            elem_type = TypeElement;
            NNode=size(CONNECT,2);
            ndime = size(COOR,2) ;
            nnod = size(COOR,1) ;
            npe =  size(CONNECT,2) ;
            nElem = size(CONNECT,1) ;
            fid = fopen(NameFile,'wt');
            fprintf(fid,' # ================================================== \n',[]);
            fprintf(fid,[' # \n'],[]);
            fprintf(fid,' # POSTPROCESSING WITH GID - MESH FILE \n',[]);
            fprintf(fid,[' # EXAMPLE NAME: ' NAMEPROJ ' \n'],[]);
            fprintf(fid,' # ================================================== \n',[]);
            fprintf(fid,[' MESH "' NAMEPROJ '" dimension ' num2str(ndime) ' Elemtype ' elem_type ...
                ' Nnode ' num2str(NNode) '\n'],[]);
            fprintf(fid,' Coordinates \n',[]);
            fprintf(fid,' # node number   coordinate_x  coordinate_y  coordinate_z \n',[]);
            format_xnod = [];
            for k=1:ndime
                format_xnod = [format_xnod, '%15.5e'];
            end
            format_xnod=[' %10i ' format_xnod, '\n'];
            fprintf(fid,format_xnod,[numer_nodes',COOR(:,1:ndime) ]');
            fprintf(fid,' end coordinates \n',[]);
            fprintf(fid,' Elements \n',[]);
            format_title = ' # element ';
            for k=1:npe
                format_title=[format_title, ' node_' num2str(k)];
            end
            format_title=[format_title, ' material number \n'];
            fprintf(fid,format_title,[]);
            format_icone = [];
            for k=1:npe
                format_icone=[format_icone, '%7i'];
            end
            format_icone=['%10i ' format_icone, '     %4i ' '\n'];
            fprintf(fid,format_icone,[(1:nElem)',CONNECT, MaterialType]');
            fprintf(fid,' end elements ',[]);
            fclose(fid);
        end

        function GidResults2DFE_modes(obj,NameFile,COOR,TypeElement,MODES,posgp)
            elem_type = TypeElement;
            ndime = size(COOR,2) ;
            nnod = size(COOR,1) ;
            npg = size(posgp,2);
            fid_res = fopen(NameFile,'wt');
            fprintf(fid_res,'GiD Post Results File 1.0 \n');
            xg =posgp ;
            fprintf(fid_res,['GaussPoints "GPset" Elemtype ',elem_type,'\n']);
            fprintf(fid_res,['Number of Gauss Points: ',num2str(npg),'\n']);
            fprintf(fid_res,'Nodes not included\n');
            fprintf(fid_res,'Natural Coordinates: Internal\n');
            fprintf(fid_res,'End GaussPoints\n');
            TIMEVECTOR = 1;
            istep = 1;
            disp(['istep=',num2str(istep)])
            for imode = 1:size(MODES,2)
                d = MODES(:,imode) ;
                var = reshape(d,ndime,nnod) ;
                time_step = TIMEVECTOR(istep) ;
                fprintf(fid_res, ['Result "Mode ',num2str(imode),'"  "Load Analysis" '  num2str(time_step) ' Vector OnNodes  \n' ],[]);
                if ndime==2
                    fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
                elseif ndime==3
                    fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
                end
                fprintf(fid_res,'Values\n');
                fORMAT = ['%d',repmat(' %f',1,ndime),'\n'];
                fprintf(fid_res,fORMAT,[1:nnod;var]);
                fprintf(fid_res,'End Values\n');
            end
            fclose(fid_res);
        end

    end

end