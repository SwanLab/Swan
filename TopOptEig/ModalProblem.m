classdef ModalProblem < handle
    
    properties (Access = public)
        variables
    end
    
    properties (Access = private)
        boundaryConditions
        freeDOFs
        stiffnessMatrix
        massMatrix
        geometry
        scale
        pdim
        ptype
        inputBC
        V
        modes
        modesF
        Mode1x
        Mode1y
    end

    properties (Access = private)
        quadrature
        mesh
        dim
        material
        dofConnec
        interpolation
        interpolationType
        displacementField
    end
    
    methods (Access = public)

        function obj = ModalProblem(cParams)
            obj.init(cParams)
            obj.computeDimensions();
            obj.createDisplacementField();
            obj.createBoundaryConditions();
            obj.createDOFsconnect();
        end

        function solve(obj,x)
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix(x);
            obj.computeEigModes();
            obj.computeStrain();
            % obj.printInGiD();
            % obj.computeStress();
        end

        function dim = getDimensions(obj)
            dim = obj.displacementField.dim;
        end

        function setC(obj, C)
            obj.material.C = C;
        end


    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.scale    = cParams.scale;
            obj.pdim     = cParams.dim;
            obj.ptype    = cParams.type;
            obj.inputBC  = cParams.bc;  
            if isprop(cParams, 'interpolationType')
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end
            obj.createQuadrature();
            obj.createInterpolation();
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,obj.interpolationType);
%             int = Interpolation.create(obj.mesh,'QUADRATIC');
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end

        function computeDimensions(obj)
            s.type      = 'Vector';
            s.fieldName = 'Disp'; %'u';
            s.mesh      = obj.mesh;
            s.ndimf = str2double(regexp(obj.pdim,'\d*','Match'));
            d = DimensionVariables.create(s);
            d.compute()
            obj.dim = d;
        end

        function createDisplacementField(obj)
            ndimf = regexp(obj.pdim,'\d*','Match');
            s.mesh               = obj.mesh;
            s.ndimf              = str2double(ndimf);
            s.inputBC            = obj.inputBC;
            s.scale              = obj.scale;
            s.interpolationOrder = obj.interpolationType; % obj.interpolationType
            obj.displacementField = Field(s);
        end

        function createBoundaryConditions(obj)
            bc = obj.displacementField.inputBC;
            bc.ndimf = obj.displacementField.dim.ndimf;
            bc.ndofs = obj.displacementField.dim.ndofs;
            s.dim   = obj.displacementField.dim;
            s.mesh  = obj.mesh;
            s.scale = 'MACRO';
            s.bc    = {bc};
            s.ndofs = obj.displacementField.dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createDOFsconnect(obj)
            connec  = obj.mesh.connec;
            ndimf   = obj.dim.ndimf;
            nnodeEl = size(connec, 2);
            ndofsEl = nnodeEl * ndimf;
            dofsElem  = zeros(ndofsEl,size(connec,1));
            for inode = 1:nnodeEl
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            obj.dofConnec = dofsElem;
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimf;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrix';
            s.mesh          = obj.mesh;
            s.globalConnec  = obj.displacementField.connec;
            s.dim           = obj.displacementField.dim;
            s.material      = obj.material;
            s.interpolation = obj.interpolation;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            obj.stiffnessMatrix = K;
        end

        function computeMassMatrix(obj,xReg)
            s.type = 'MassMatrix';
            s.mesh          = obj.mesh;
            s.quadType      = 'QUADRATICMASS';
            s.globalConnec  = obj.displacementField.connec;
            s.dim           = obj.displacementField.dim;
            s.material      = obj.material;
            s.interpolation = obj.interpolation;
            s.quadrature    = obj.quadrature;            
            LHS = LHSintegrator.create(s);
            M   = LHS.compute(); % xReg 
            obj.massMatrix = M;
        end

        function computeEigModes(obj)
            bc = obj.boundaryConditions;
            K = obj.stiffnessMatrix;
            M = obj.massMatrix;
            Kfree = bc.fullToReducedMatrix(K);
            Mfree = bc.fullToReducedMatrix(M);
            [v,d] = eigs(Kfree,Mfree,2,'SM');
            obj.computeBucklingModes(v);
            obj.variables.eigenModes  = obj.V;
            lambda = obj.computeLambda(d);
            obj.variables.eigenValues  = lambda;
        end

        function computeBucklingModes(obj,v)
            obj.modesF = v
            ndof = obj.dim.ndofs;
            Modes=zeros(ndof,2);
            free = obj.boundaryConditions.free;
            obj.freeDOFs = free;
            Modes(free,1) = v(:,1);
            Modes(free,2) = v(:,2);
            obj.Mode1x = Modes(1:2:end-1,1);
            obj.Mode1y = Modes(2:2:end,1);
            obj.variables.modes = Modes;
            obj.V = Modes;
        end

        function l = computeLambda(obj,d)
            l = sort(diag(d));
        end

        function computeStrain(obj)
            %obj.createDisplacementField();
            s.dim          = obj.displacementField.dim;
            s.mesh         = obj.mesh;
            s.quadrature   = obj.quadrature;
            s.displacement = obj.variables.eigenModes(:,1);
            s.dispField    = obj.displacementField;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;        
        end

        function printInGiD(obj)
            fileName = 'EigModesProdois';
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
            % f.u = obj.mode1;
            dP.fields.u = [obj.Mode1x obj.Mode1y];
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

%         function printInGiD(obj)
%             coor = obj.mesh.coord;
%             conn = obj.mesh.connec;
%             MODES = obj.variables.eigenModes;
%             posgp = 1;
%             DOFl = obj.freeDOFs;
%             TypeElement = 'TRIANGLE';
%             NameFileMesh = 'PostModes';
%             NameFile_msh = ['GIDPOST/','MODES2','_',NameFileMesh,'.msh'] ;
%             NameFile_res= ['GIDPOST/','MODES2','_',NameFileMesh,'.res'] ;
%             MaterialType = ones(size(conn,1),1) ; 
%             obj.GidMesh2DFE(NameFile_msh,coor,conn,MaterialType,TypeElement);
%             MODESplot = zeros(length(obj.dim.ndofs),size(MODES,2)) ;
%             MODESplot = obj.V;% (DOFl,:) 
%             obj.GidResults2DFE_modes(NameFile_res,coor,TypeElement,MODESplot,posgp);
%             cddd = cd ;
%             NAMEFILEOPEN =  [cddd,'/',NameFile_res] ;
%             disp('open GID FILE:')
%             disp(NAMEFILEOPEN)
%         end

         function GidMesh2DFE(obj,NameFile,COOR,CONNECT,MaterialType,TypeElement)
            NAMEPROJ = 'MODES2';
            numer_nodes = 1:size(COOR,1);
            elem_type = TypeElement;
            NNode=size(CONNECT,2);
            ndime = size(COOR,2) ;
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