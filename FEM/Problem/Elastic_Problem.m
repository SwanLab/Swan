classdef Elastic_Problem < FEM
    %Elastic_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    properties (Access = private)
        material
        fileName
        nFields
        interp
    end
    
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem(fileName)
            obj.fileName = fileName;
            obj.nFields = 1;
            obj.readProblemData(fileName);
            obj.createGeometry();
            obj.createInterpolation();
            obj.createDOF();
            obj.createMaterial();
            obj.createSolver(); 
            obj.createElement();            
        end

        function computeVariables(obj)
            Kred = obj.element.computeLHS();
            obj.element.computeRHS();
            u = obj.solver.solve(Kred,obj.element.fextRed);
            obj.variables = obj.element.computeVars(u);
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

    end
    
    methods (Access = private)
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
        end     
        
        function createInterpolation(obj)
            obj.interp{1} = Interpolation.create(obj.mesh,'LINEAR');                                    
        end        
        
        function createMaterial(obj)
            s.ptype = obj.problemData.ptype;
            s.pdim  = obj.problemData.pdim;
            s.nelem = obj.mesh.nelem;
            s.geometry = obj.geometry;
            s.mesh  = obj.mesh;            
            obj.material = Material.create(s);                        
        end
        
        function createSolver(obj)
            obj.solver = Solver.create;
        end
        
        function createDOF(obj)
            obj.dof = DOF_Elastic(obj.fileName,obj.mesh,obj.problemData.pdim,obj.nFields,obj.interp);            
        end
        
        function createElement(obj)
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,obj.material,obj.dof,obj.problemData,obj.interp);            
        end
        
    end
    
end

