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
        interp
    end
    
    % Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem_Micro(fileName)
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
        
        function setMesh(obj,mesh)
            obj.mesh = mesh;
            obj.createGeometry();
            obj.createInterpolation();            
            obj.setDOF();
            obj.createMaterial();
            obj.createSolver();
            obj.createElement();            
        end
        
        function [Ch,tstrain,tstress] = computeChomog(obj)
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
            nstre = obj.element.getNstre();
            basis = diag(ones(nstre,1));
            tstrain  = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tstrainF = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tstress  = zeros(nstre,obj.element.quadrature.ngaus,nstre,obj.element.nelem);
            tdisp   = zeros(nstre,obj.element.dof.ndof);
            
            Ch2 = zeros(nstre,nstre);
            Ch = zeros(nstre,nstre);
            
            var2print = cell(nstre,1);
            for istre=1:nstre
                strain = basis(istre,:);
                obj.element.setVstrain(strain);
                obj.computeVariables();
                Ch(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstrainF(istre,:,:,:) = obj.variables.strain_fluct;
                tstress(istre,:,:,:) = obj.variables.stress;
                
                tdisp(istre,:) = obj.variables.d_u;
                var2print{istre}.stress = obj.variables.stress;
                var2print{istre}.strain = obj.variables.strain;
                var2print{istre}.stress_fluct = obj.variables.strain_fluct;
                var2print{istre}.strain_fluct = obj.variables.strain_fluct;
                var2print{istre}.d_u = obj.variables.d_u;
                var2print{istre}.fext = obj.variables.fext;
            end
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
            obj.variables.tdisp   = tdisp;
            
            dV = obj.geometry.dvolu;
            for istre = 1:nstre
                for jstre = 1:nstre
                    s = squeezeParticular(tstress(istre,:,:,:),1);
                    e = squeezeParticular(tstrain(jstre,:,:,:),1);
                    ener = (s.*e);
                    en = sum(ener,2);
                    en = squeezeParticular(en,2);
                    Ch2(istre,jstre) = en(:)'*dV(:);
                end
            end
            
            
            Cmat = obj.element.material.C;
            
            Ch3 = zeros(nstre,nstre);
            ngaus = size(tstrain,2);
            nelem = size(tstrain,4);
            for istre = 1:nstre
                for jstre = 1:nstre
                   ei(1:ngaus,:,:) = squeeze(tstrain(istre,:,:,:));
                   ej(1:ngaus,:,:) = squeeze(tstrain(jstre,:,:,:));
                   c = zeros(ngaus,nelem);
                    for kstre = 1:nstre
                        for lstre = 1:nstre
                         eiV(1:ngaus,:) = squeeze(ei(:,kstre,:));
                         ejV(1:ngaus,:) = squeeze(ej(:,lstre,:));
                         Cm = squeezeParticular(Cmat(kstre,lstre,:,:),1);
                         Cm = squeezeParticular(Cm,1)';
                         %Cm  = repmat(Cmm,ngaus,1);
                         c = c + Cm.*eiV.*ejV;
                        end
                    end
                    cC = c.*dV';                    
                    Ch3(istre,jstre) = sum(cC(:));
                end
            end
            

            
            obj.variables2print = var2print;
            
            obj.Chomog = Ch;
            obj.tstrain = tstrain;
            obj.tstress = tstress;
        end
        
        function computeStressBasisCellProblem(obj)
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
            nstre = obj.element.getNstre();
            basis = diag(ones(nstre,1));          
            var2print = cell(nstre,1);
            for istre=1:nstre
                stress = basis(istre,:);
                v = obj.computeVariablesForGivenStress(stress,obj.Chomog);
                var2print{istre}= v;
            end
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
        
        function stress = computeVarFromStress(obj,stress)
           Ch = obj.computeChomog();
           varFromCh = obj.variables;
           v  = obj.computeVariablesForGivenStress(stress,Ch);  
           v.varFromCh = varFromCh;
           obj.variables = v;
           stress = v.stress;
        end
        
        function sPnorm = computeStressPnorm(obj,stress,p)
            stresses = obj.computeVarFromStress(stress);
            m = obj.mesh;
            quad = Quadrature.set(m.geometryType);
            quad.computeQuadrature('CONSTANT');
            dV = m.computeDvolume(quad);
            sx  = squeeze(stresses(:,1,:));
            sy  = squeeze(stresses(:,2,:));
            sxy = squeeze(stresses(:,3,:));
            sNorm = sqrt(sx.*sx + 2*sxy.*sxy + sy.*sy);
            if isequal(p,'max')
                sPnorm = max(sNorm);
            else
                int = sNorm.^p;
                sPnorm = sum(int(:).*dV(:))^(1/p);
            end            
        end
        
    end
    
    methods (Access = private)
        
        function v = computeVariablesForGivenStress(obj,stress,Chomog)
            strain(1,:) = Chomog\stress';
            obj.element.setVstrain(strain);
            obj.computeVariables();
            v.stress = obj.variables.stress;
            v.strain = obj.variables.strain;
            v.stress_fluct = obj.variables.strain_fluct;
            v.strain_fluct = obj.variables.strain_fluct;
            v.d_u = obj.variables.d_u;
            v.fext = obj.variables.fext;
            v.Ch   = obj.Chomog;
        end        
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
        end
        
        function createInterpolation(obj)
            obj.interp{1} = Interpolation.create(obj.mesh,'LINEAR');                                    
        end
        
        function createMaterial(obj)
            cParams.ptype = obj.problemData.ptype;
            cParams.pdim  = obj.problemData.pdim;
            cParams.nelem = obj.mesh.nelem;
            cParams.geometry = obj.geometry;
            cParams.mesh  = obj.mesh;
            obj.material = Material.create(cParams);
        end
        
        function createSolver(obj)
            obj.solver = Solver.create;
        end
        
        function createDOF(obj)
            obj.dof = DOF_Elastic_Micro(obj.fileName,obj.mesh,obj.problemData.pdim,obj.nFields,obj.interp);
        end
        
        function setDOF(obj,mesh)
            obj.dof.setMesh(obj.mesh);
        end        
        
        function createElement(obj)
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,obj.material,obj.dof,obj.problemData,obj.interp);
        end
        
    end
end