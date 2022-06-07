classdef ModalProblem < handle
    
    properties (Access = public)
        variables
        dim
        mesh
        dofConnec
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
    end
    
    properties (Access = private)
        quadrature
        
        material

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
            % obj.createSolver();
        end

        function solve(obj,x)
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix(x);
            obj.computeEigModes();
            obj.computeStrain();
            obj.printInGiD();
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
            obj.inputBC  = cParams.bc;   % No necesario, se calculan aquí directo
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
            s.fieldName = 'Disp';%'u';
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

        %         function createBoundaryConditions(obj)
        %             d = obj.dim;
        %             FixNod = obj.computeFixedNodes();
%             FixDof = obj.computeFixedDOFs(FixNod); 
%             dofs = 1:d.ndofs;
%             free  = setdiff(dofs,FixDof);
%             obj.freeDOFs = free';
%         end  

%         function FN = computeFixedNodes(obj)
%             coorX = obj.mesh.coord(:,1);
%             % coorY = obj.mesh.coord(:,2);
%             tol = 0.001;
%             minX = min(coorX) + tol;
%             maxX = max(coorX) - tol;
%             FN1 = find(coorX < minX);
%             FN2 = find(coorX > maxX);
%             FN = [FN1 ; FN2];
%         end

%         function FixDof = computeFixedDOFs(obj, FixNod)
%             d = obj.dim;
%             nDOFn = d.ndofs/d.nnodes;
%             FD = zeros(nDOFn*length(FixNod),1);
%             for i = 1: length(FixNod)
%                 FD(2*i-1,1) = FixNod(i)*2-1;
%                 FD(2*i,1) = FixNod(i)*2;
%             end
%             FixDof = FD';
%         end

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
            s.type = 'MassMatrixModal';
            s.mesh          = obj.mesh;
            s.quadType      = 'QUADRATICMASS';
            s.globalConnec  = obj.displacementField.connec;
            s.dim           = obj.displacementField.dim;
            s.material      = obj.material;
            s.interpolation = obj.interpolation;
            s.quadrature    = obj.quadrature;            
            LHS = LHSintegrator.create(s);
            M   = LHS.compute(xReg); % xReg 
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
            %obj.variables.mode1 = obj.V(:,1);
            %obj.variables.mode2 = obj.V(:,2);
            lambda = obj.computeLambda(d);
            obj.variables.eigenValues  = lambda;
        end

        function computeBucklingModes(obj,v)
            ndof = obj.dim.ndofs;
            Modes=zeros(ndof,2);
            free = obj.boundaryConditions.free;
            Modes(free,1) = v(:,1);
            Modes(free,2) = v(:,2);
            obj.variables.modes = Modes;
            obj.V = Modes;
        end

%         function MatrixFree = provideFreeMatrix(obj,Matrix)
%             free = obj.freeDOFs;
%             MatrixFree = Matrix(free,free);
%         end

%         function V = addBoundaryConditions(obj,v)
%             free = obj.freeDOFs;
%             obj.V = zeros(obj.dim.ndofs,2);
%             obj.V(free,1) = v(:,1);
%             obj.V(free,2) = v(:,2);
%             V = obj.V;
%         end

%         function V_f = filterDOFtoELEM(obj,V)
%             connec = obj.mesh.connec;
%             nnodeEl = obj.dim.nnodeElem;
%             nunkn   = obj.dim.ndimf;
%             for inode=1:nnodeEl
%                 nodes = connec(:,inode);
%                 for idime = 1:nunkn
%                     dofs = nunkn*(nodes - 1) + idime;
%                     V_f = V(dofs,:);
%                 end
%             end
%         end

        function l = computeLambda(obj,d)
            l = sort(diag(d));
        end

        function computeStrain(obj)
            obj.createDisplacementField();
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
            fileName = 'EigModesPro2';
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
            dP.fields.u = obj.variables.modes;
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