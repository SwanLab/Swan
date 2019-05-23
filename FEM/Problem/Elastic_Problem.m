classdef Elastic_Problem < FEM
    %Elastic_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Elastic(problemID,obj.geometry,obj.mesh);
        end
        
        function preProcess(obj)
            cParams.ptype = obj.mesh.ptype;
            cParams.pdim  = obj.mesh.pdim;
            cParams.nelem = obj.geometry(1).interpolation.nelem;
            cParams.geometry = obj.geometry;
            cParams.mesh  = obj.mesh;            
            material = Material.create(cParams);
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,material,obj.dof);
            obj.solver = Solver.create;
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
            f = obj.element.bcApplier.full_vector_2_reduced_vector(fbody);
            u = obj.solver.solve(Kred,f);
            obj.variables = obj.element.computeVars(u);                        
        end
        
        function c = computeCompliance(obj)
            dvolum  = obj.geometry.dvolu;
            strain = obj.variables.strain;
            C = obj.element.material.C;
            ngaus = obj.element.quadrature.ngaus;
            nstre = size(strain,2);
            c = 0;
            for igaus = 1:ngaus
                strainG = squeeze(strain(igaus,:,:));
                dV = dvolum(:,igaus);
                for istre = 1:nstre
                    sI(:,1) = strainG(istre,:);
                    for jstre = 1:nstre
                        sJ(:,1) = strainG(jstre,:);
                        Cij = squeeze(C(istre,jstre,:));
                        csum = sI.*Cij.*sJ.*dV;
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
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
    end
end

