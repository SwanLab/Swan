classdef Element_Stokes < Element
    %Element_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS_elem
        LHS
        RHS
    end

    properties(Access = private)
        D_elem
        M_elem
        K_elem
        D
        dt
        mesh
        velocityField
        pressureField
        forcesFormula
        LHSint
        boundaryConditions
    end
    
    methods
        function obj = Element_Stokes(geom,mesh,material,dof,pD,interp, ...
                vField, pField,forcesFormula, BC)
            obj.initElement(geom,mesh,material,dof,pD.scale,interp);
%             obj.dim = dim;
            obj.mesh = mesh;
            %obj.nstre=0;
            obj.nfields=2;
            obj.velocityField = vField;
            obj.pressureField = pField;
            obj.forcesFormula = forcesFormula;
            obj.boundaryConditions = BC;
        end
        
        function [r,dr] = computeResidual(obj,x,dr,x_n)
%             K = compute_LHS(obj);
            freeV = obj.boundaryConditions.freeFields{1};
            lenFreeV = length(freeV);
            if (nargin ==3)
                % Steady
                Mred_x_n = zeros(lenFreeV,1);
            else
                % Transient
                Melem = obj.LHSint.Melem;
                s.dim          = obj.velocityField.dim;
                s.globalConnec = obj.velocityField.connec;
                s.nnodeEl      = obj.velocityField.dim.nnodeElem;
                assembler = Assembler(s);
                M = assembler.assemble(Melem);
                M = obj.symGradient(M);

                Mred = M(freeV,freeV);
                Mred_x_n = Mred*x_n;
            end
            
            Fext = compute_RHS(obj);
            
            
            R = obj.compute_imposed_displacement_force(obj.LHS);
            Fext = Fext + R ;
            
            
%             Fext_red = obj.bcApplier.fullToReducedVector(Fext);
            Fext_red = obj.boundaryConditions.fullToReducedVector(Fext);
            Fext_red(1:lenFreeV,1) = Fext_red(1:lenFreeV,1) + Mred_x_n;
            
            fint_red = dr*x;
            
            r = fint_red - (Fext_red);
%             dr = Kred;
            
        end
        
        function R = compute_imposed_displacement_force(obj,K)
            % Forces coming from imposed displacement
            dirichlet = obj.boundaryConditions.dirichlet;
            uD = obj.boundaryConditions.dirichlet_values;
            R  = -K(:,dirichlet)*uD;
        end
        
        function dr = computedr(obj,dt)
            if nargin < 2
                dt=inf;
            end
            obj.LHS = compute_LHS(obj,dt);
%             LHSred = obj.bcApplier.fullToReducedMatrix(obj.LHS);
            LHSred = obj.boundaryConditions.fullToReducedMatrix(obj.LHS);
            dr = LHSred;
        end
        
        function LHS = compute_LHS(obj,dt)
%             obj.dt = dt;
%             AA = obj.computeVelocityLaplacian();
%             D = obj.computeDmatrix();
%             BB = obj.computePressureLHSMatrix();
%             LHS = [AA, D; D',BB];

            % Inefficient. It is always the same.
            s.type          = 'Stokes';
            s.dt            = dt;
            s.mesh          = obj.mesh;
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            LHS_int = LHSintegrator.create(s);
            LHS = LHS_int.compute();
            obj.LHSint = LHS_int;
        end
        
        function RHS = compute_RHS(obj)
            % Inefficient. It is always the same.
            s.type          = 'Stokes';
            s.mesh          = obj.mesh;
            s.velocityField = obj.velocityField;
            s.pressureField = obj.pressureField;
            s.forcesFormula = obj.forcesFormula;
            RHS_int = RHSintegrator.create(s);
            RHS_elem = RHS_int.integrate();
            RHS = AssembleVector(obj,RHS_elem);
        end
        
        function variable = computeVars(obj,x_free)
%             x = obj.bcApplier.reducedToFullVector(x_free);
            x = obj.boundaryConditions.reducedToFullVector(x_free);
            ndofsV = obj.velocityField.dim.ndofs;
            variable.u = x(1:ndofsV,:);
            variable.p = x(ndofsV+1:end,:);
        end
    end

    methods (Access = private)

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end
    end

end