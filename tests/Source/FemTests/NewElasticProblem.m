classdef NewElasticProblem < NewFEM
    %Elastic_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        meshNova
        lhs
        integrator
    end

    properties (Access = private)
        material
        fileName
        nFields
        interp
        bcApplier
    end

    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = NewElasticProblem(fileName)
            obj.fileName = fileName;
            obj.nFields = 1;
            obj.readProblemData(fileName); % creem Mesh
            %%
            % Ara amb integrators
            obj.NewCreateIntegrators();
            %%
%             obj.createGeometry(); % inclou quadratura (lhsint) Autoconsum per DOF i element
%             obj.createInterpolation(); %lhsintegrator
%             obj.createDOF();
%             obj.createMaterial();
            obj.createSolver(); 
%             obj.createElement();  % no te counterpart  
        end

        function [K, Kred] = getElementLHS(obj)
            K = obj.element.K;
            Kred = obj.element.computeLHS();
        end

        function el = getElement(obj)
            el = obj.element;
        end

        function computeVariables(obj)
%             Kred = obj.element.computeLHS();
%             obj.element.computeRHS();
%             u = obj.solver.solve(Kred,obj.element.fextRed);
%             obj.variables = obj.element.computeVars(u);

            Kred = obj.integrator.Kred;
            fNodal = obj.computeExternalForces();
%             R = obj.compute_imposed_displacement_force(obj.K);
%             obj.fext = Fext + R;
            fext_red = obj.bcApplier.fullToReducedVector(fNodal);
%             obj.rhs = obj.integrator.integrate(fNodal);
            u = obj.solver.solve(Kred,fext_red);
            obj.variables = obj.computeVars(u);
        end

        function computeVariablesWithBodyForces(obj,fbody)
            Kred = obj.element.computeLHS();
            %obj.element.computeRHS();
            f = obj.element.bcApplier.fullToReducedVector(fbody);
            u = obj.solver.solve(Kred,f);
            obj.variables = obj.element.computeVars(u);
        end
        
        function c = computeCompliance(obj)
            dvolum  = obj.geometry.dvolu;
            stress = obj.variables.stress;
            C = obj.element.material.C;
            
            C11 = C(1,1,:);
            C22 = C(2,2,:);
            C12 = C(1,2,:);
            C33 = C(3,3,:);
            
            Cinv = zeros(size(C));

            det = C11.*C22 - C12.^2;
            Cinv(1,1,:) = C22./det;
            Cinv(2,2,:) = C11./det;
            Cinv(1,2,:) = -C12./det;
            Cinv(2,1,:) = -C12./det;
            Cinv(3,3,:) = 1./C33;

            ngaus = obj.element.quadrature.ngaus;
            nstre = size(stress,2);
            c = 0;
            for igaus = 1:ngaus
                stressG = squeeze(stress(igaus,:,:));
                dV = dvolum(:,igaus);
                for istre = 1:nstre
                    sI(:,1) = stressG(istre,:);
                    for jstre = 1:nstre
                        sJ(:,1) = stressG(jstre,:);
                        %Cij = squeeze(C(istre,jstre,:));
                        %csum = sI.*Cij.*sJ.*dV;
                        
                        Cinv_ij = squeeze(Cinv(istre,jstre,:));
                        csum = sI.*Cinv_ij.*sJ.*dV;
                        
                        c = c + sum(csum);
                    end
                end
            end
        end
        
%         function print(obj)
%             postprocess = Postprocess_PhysicalProblem;
%             results.physicalVars = obj.variables;
%             postprocess.print(obj,obj.problemID,results);
%         end
        
        function postProcess(obj)
            % ToDo
            % Inspire in TopOpt
        end
        
        % !! THIS SHOULD BE DEFINED BY THE USER !!
        function NewCreateIntegrators(obj)
            s.type = 'SIMPLE';
            s.mesh  = obj.mesh;
            s.npnod = obj.mesh.npnod;
            s.fileName     = obj.fileName;
            s.globalConnec = obj.mesh.connec;
            s.problemData  = obj.problemData;
            obj.integrator = Integrator.create(s);
            obj.lhs = obj.integrator.computeLHS();
            obj.bcApplier = obj.integrator.bcApplier;
            obj.dof = obj.integrator.dof;
            obj.interp{1} = Interpolation.create(obj.mesh,'LINEAR');
        end

    end
    
    methods (Access = private)

%         function createGeometry(obj)
%             s.mesh = obj.mesh;
%             obj.geometry = Geometry.create(s);
%         end     
%         
%         function createInterpolation(obj)
%             obj.interp{1} = Interpolation.create(obj.mesh,'LINEAR');
%         end        
%         
%         function createMaterial(obj)
%             s.ptype = obj.problemData.ptype;
%             s.pdim  = obj.problemData.pdim;
%             s.nelem = obj.mesh.nelem;
%             s.geometry = obj.geometry;
%             s.mesh  = obj.mesh;            
%             obj.material = Material.create(s);
%         end
        
        function createSolver(obj)
            obj.solver = Solver.create;
        end
        
%         function createDOF(obj)
%             obj.dof = DOF_Elastic(obj.fileName,obj.mesh,obj.problemData.pdim,obj.nFields,obj.interp);            
%         end
%         
%         function createElement(obj)
%             obj.element = Element_Elastic.create(obj.mesh,obj.geometry,obj.material,obj.dof,obj.problemData,obj.interp);            
%         end
        
        
        function variables = computeVars(obj, uL)
            variables.d_u = obj.computeDisplacements(uL);
%             variables.fext = obj.fext;
%             variables.strain = obj.computeStrain(variables.d_u,obj.dof.in_elem{1});
%             variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.quadrature.ngaus,obj.nstre);
%             variables = obj.permuteStressStrain(variables);
%             if ~isequal(class(obj),'Element_Elastic_3D')
%                 [dir,s] = obj.computePrincipalStressDirection(variables.stress);
%                 variables.principalDirections = dir;
%                 variables.principalStress     = s;
%             end
        end

        function u = computeDisplacements(obj, usol)
            u = obj.bcApplier.reducedToFullVector(usol);
        end

        function Fext = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext;
            FextVolumetric  = obj.computeVolumetricFext;
            FextSupVol = {FextSuperficial + FextVolumetric};
            FextSupVol = obj.AssembleVector(FextSupVol);
            FextPoint = obj.computePunctualFext();
            Fext = FextSupVol +  FextPoint;
        end

        function FextSuperficial = computeSuperficialFext(obj)
            nnode = 3;
            nunkn = 2;
            nelem = 16;
            FextSuperficial = zeros(nnode*nunkn,1,nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            nnode = 3;
            nunkn = 2;
            nelem = 16;
            FextVolumetric = zeros(nnode*nunkn,1,nelem);
        end

        function b = AssembleVector(obj,b_elem_cell)
            nfields = 1;
            for ifield = 1:nfields
                b_elem = b_elem_cell{ifield,1};
                b = zeros(obj.dof.ndof(ifield),1);
                for i = 1:obj.interp{ifield}.nnode*obj.dof.nunkn(ifield)
                    for igaus = 1:size(b_elem,2)
                    c = squeeze(b_elem(i,igaus,:));
                    idof_elem = obj.dof.in_elem{ifield}(i,:);
                    b = b + sparse(idof_elem,1,c',obj.dof.ndof(ifield),1);
                    end
                end
                b_global{ifield,1} = b;
            end
            b=cell2mat(b_global);
        end

        function FextPoint = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            FextPoint = zeros(obj.dof.ndof,1);
            if ~isempty(obj.dof.neumann)
                FextPoint(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end
    end

end