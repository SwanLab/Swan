classdef Elastic_Problem_Micro < FEM
    %Elastic_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    % Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        variables2printStressBasis
        variables2print
    end
    
    % Private properties definition ======================================
    properties (Access = private)
        Chomog
        tstrain
        tstress
        
        material
        fileName
        nFields
    end
    
    % Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem_Micro(fileName)
            obj.fileName = fileName;
            obj.nFields = 1;
            obj.readProblemData(fileName);
            obj.createGeometry();
            obj.createDOF();
            obj.createMaterial();
            obj.createSolver();
            obj.createElement();
        end
        
        function [Chomog,tstrain,tstress] = computeChomog(obj)
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.interpolation_u.computeShapeDeriv(obj.element.quadrature.posgp)
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
            nstre = obj.element.getNstre();
            basis = diag(ones(nstre,1));
            tstrain = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tstress = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tdisp   = zeros(nstre,obj.element.dof.ndof);
            
            var2print = cell(nstre,1);
            for istre=1:nstre
                strain = basis(istre,:);
                obj.element.setVstrain(strain);
                obj.computeVariables;
                Chomog(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstress(istre,:,:,:) = obj.variables.stress;
                tdisp(istre,:) = obj.variables.d_u;
                var2print{istre}.stress = obj.variables.stress;
                var2print{istre}.strain = obj.variables.strain;
                var2print{istre}.stress_fluct = obj.variables.strain_fluct;
                var2print{istre}.strain_fluct = obj.variables.strain_fluct;
                var2print{istre}.d_u = obj.variables.d_u;
                var2print{istre}.fext = obj.variables.fext;
            end
            obj.variables.Chomog  = Chomog;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
            obj.variables.tdisp   = tdisp;
            
            
            obj.variables2print = var2print;
            
            obj.Chomog = Chomog;
            obj.tstrain = tstrain;
            obj.tstress = tstress;
        end
        
        function computeStressBasisCellProblem(obj)
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.interpolation_u.computeShapeDeriv(obj.element.quadrature.posgp)
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
            nstre = obj.element.getNstre();
            basis = diag(ones(nstre,1));
            Chomog =  zeros(nstre,nstre);
            tstrain = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tstress = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tdisp   = zeros(nstre,obj.element.dof.ndof);
            var2print = cell(nstre,1);
            for istre=1:nstre
                stress = basis(istre,:);
                strain(1,:) = (obj.Chomog)\stress';
                obj.element.setVstrain(strain);
                obj.computeVariables;
                Chomog(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstress(istre,:,:,:) = obj.variables.stress;
                tdisp(istre,:) = obj.variables.d_u;
                var2print{istre}.stress = obj.variables.stress;
                var2print{istre}.strain = obj.variables.strain;
                var2print{istre}.stress_fluct = obj.variables.strain_fluct;
                var2print{istre}.strain_fluct = obj.variables.strain_fluct;
                var2print{istre}.d_u = obj.variables.d_u;
                var2print{istre}.fext = obj.variables.fext;
            end
            %            obj.variables2printStressBasis.Chomog  = Chomog;
            % obj.variables2printStressBasis.tstrain = tstrain;
            %  obj.variables2printStressBasis.tstress = tstress;
            %   obj.variables2printStressBasis.tdisp   = tdisp;
            obj.variables2printStressBasis = var2print;
            
        end
        
        
        function v = computeGeometricalVolume(obj)
            v = 1;%sum(sum(obj.geometry.dvolu));
        end
        
        function computeVariables(obj)
            Kred = obj.element.computeLHS();
            obj.element.computeRHS();
            u = obj.solver.solve(Kred,obj.element.fextRed);
            obj.variables = obj.element.computeVars(u);
        end        
        
        
        
    end
    
    methods (Access = private)
        
        function createGeometry(obj)
            obj.geometry = Geometry(obj.mesh,'LINEAR');
        end
        
        function createMaterial(obj)
            cParams.ptype = obj.problemData.ptype;
            cParams.pdim  = obj.problemData.pdim;
            cParams.nelem = obj.geometry(1).interpolation.nelem;
            cParams.geometry = obj.geometry;
            cParams.mesh  = obj.mesh;
            obj.material = Material.create(cParams);
        end
        
        function createSolver(obj)
            obj.solver = Solver.create;
        end
        
        function createDOF(obj)
            obj.dof = DOF_Elastic_Micro(obj.fileName,obj.geometry,obj.problemData.pdim,obj.nFields);
        end
        
        function createElement(obj)
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,obj.material,obj.dof,obj.problemData);
        end
        
    end
end