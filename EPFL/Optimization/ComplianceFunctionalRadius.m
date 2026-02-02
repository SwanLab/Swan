classdef ComplianceFunctionalRadius < handle

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        stateProblem
        value0
    end

    methods (Access = public)
        function obj = ComplianceFunctionalRadius(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,mu)
%             [u,uC] = obj.computeStateVariable(mu.fun.fValues);
            [u,uC] =   obj.computeStateVariable(mu.fun);
            dK     = obj.stateProblem.computeGradK(mu.fun);
            J      = obj.computeFunction(uC);
            dJ     = obj.computeGradient(dK,uC,mu);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end


    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.stateProblem = cParams.stateProblem;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function [u,uC] = computeStateVariable(obj,mu)
            obj.stateProblem.computeLHS(mu);
%             obj.stateProblem.updateDownscaling(mu)
%             [u,uC]= obj.stateProblem.solve();
            uC = obj.stateProblem.coarseSolve();
            u=0;
            %             u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj,u)
            %             dCompliance = ElasticEnergyDensity(C,u);
            %             J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);

            % I do u'*F because i don't need to create the mesh this way.
            % It is true that to get uFine there is a reconstruction, but since
            % we are only interested in the nodal values, can use the
            % orignial mesh to get the vector.
            %continuous reconstructed
%             J = u'*obj.stateProblem.Fext;
            %discontinuous reconstructed
%              J = u(:)'*obj.stateProblem.FextDisc(:);
             %coarse
            J = u'*obj.stateProblem.EIFEMsolver.Fcoarse;
        end

        function dJ = computeGradient(obj,dK,u,mu)
            uL = obj.stateProblem.EIFEMsolver.global2local(u);
%             nelem = size(uL,2);
            [nDof, nElem] = size(uL);
            u_pages = reshape(uL, nDof, 1, nElem);
            Ku = pagemtimes(dK, u_pages);
            dj = -reshape(sum(u_pages .* Ku, 1), nElem, 1);
%             for ielem = 1:nelem
%                 dj(ielem,1) = -uL(:,ielem)'*dK(:,:,ielem)*uL(:,ielem);
%             end
            s.mesh = mu.fun.mesh;
            s.order = mu.fun.order;
            s.fValues = dj;
            dJ = {LagrangianFunction(s)};
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;
            for i = 1:length(dx)
%                 dx{i}.setFValues(dx{i}.fValues/(refX));
%                 dx{i}.setFValues(dx{i}.fValues);
                dx{i}.setFValues(dx{i}.fValues/norm(dx{i}.fValues));
            end
        end

    end

    methods (Static, Access = private)
        %         function dj = computeGradient(dK,u)
        %             uL = obj.stateProblem.EIFEMsolver.global2local(u);
        %             nelem = size(uL,2);
        %             for ielem = 1:nelem
        %                 dj(ielem) = uL(:,ielem)'*dK(:,:,ielem)*uL(:,ielem)
        %             end
        %             dj = u'*dK*u;
        %         end

       
    end

     methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end