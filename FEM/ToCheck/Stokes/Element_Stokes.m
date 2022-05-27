classdef Element_Stokes < handle
    %Element_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS
        RHS
    end

    properties(Access = private)
        dt
        mesh
        velocityField
        pressureField
        forcesFormula
        LHSint
        boundaryConditions
        material
    end
    
    methods
        function obj = Element_Stokes(mesh, material, ...
                vField, pField,forcesFormula, BC)
            obj.mesh = mesh;
            obj.material      = material;
            obj.velocityField = vField;
            obj.pressureField = pField;
            obj.forcesFormula = forcesFormula;
            obj.boundaryConditions = BC;
            obj.compute_LHS(0.01);
        end
        
        function r = computeResidual(obj,x,dr,x_n)
            freeV = obj.boundaryConditions.freeFields{1};
            lenFreeV = length(freeV);
            if (nargin ==3)
                % Steady
                Mred_x_n = zeros(lenFreeV,1);
            else
                % Transient
                % function "updateRHS"
                M = obj.LHSint.M;
                Mred = M(freeV,freeV);
                Mred_x_n = Mred*x_n;
            end
            
            F = compute_RHS(obj);
            
            
            R = obj.compute_imposed_displacement_force(obj.LHS);
            Fext = F + R ;
            
            Fext_red = obj.boundaryConditions.fullToReducedVector(Fext);
            Fext_red(1:lenFreeV,1) = Fext_red(1:lenFreeV,1) + Mred_x_n;
            
            fint_red = dr*x;
            
            r = fint_red - (Fext_red);
        end
        
        function dr = computedr(obj,dt)
            obj.LHS = compute_LHS(obj,dt);
            LHSred = obj.boundaryConditions.fullToReducedMatrix(obj.LHS);
            dr = LHSred;
        end
        
        function variable = computeVars(obj,x_free)
            x = obj.boundaryConditions.reducedToFullVector(x_free);
            ndofsV = obj.velocityField.dim.ndofs;
            variable.u = x(1:ndofsV,:);
            variable.p = x(ndofsV+1:end,:);
        end
    end

    methods (Access = private)
        
        function R = compute_imposed_displacement_force(obj,K)
            % Forces coming from imposed displacement
            dirichlet = obj.boundaryConditions.dirichlet;
            uD = obj.boundaryConditions.dirichlet_values;
            R  = -K(:,dirichlet)*uD;
        end
        
        function LHS = compute_LHS(obj,dt)
            % Inefficient. It is always the same.
            s.type          = 'Stokes';
            s.dt            = dt;
            s.mesh          = obj.mesh;
            s.material      = obj.material;
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            LHS_int = LHSintegrator.create(s);
            LHS = LHS_int.compute();
            obj.LHSint = LHS_int;
            obj.LHS = LHS;
        end
        
        function RHS = compute_RHS(obj)
            % Inefficient. It is always the same.
            s.type          = 'Stokes';
            s.mesh          = obj.mesh;
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            s.forcesFormula = obj.forcesFormula;
            RHSint = RHSintegrator.create(s);
            RHS = RHSint.integrate();
        end

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end
    end

end