classdef EigModesPlotter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        iter
        E1
        E2
        mode1
        mode2
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = EigModesPlotter(cParams)
            obj.init(cParams)            
        end
        
        function plot(obj,A,m1,m2,iter,D)
           obj.storeConvergenceVariables(iter,D)
           obj.plotColumnArea(A);
           obj.plotBucklingModes(m1,m2);
           obj.plotEigenvaluesAndIterations();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
        end

        function storeConvergenceVariables(obj,iter,D)
            obj.iter     = iter;
            obj.E1(iter) = D(1,1);
            obj.E2(iter) = D(2,2);
        end
        
        function plotColumnArea(obj,A)      
            % 1. habría que pasarle la clase dim
            z = sqrt(A); 
            nelem = obj.mesh.nelem;
            coord = obj.mesh.coord;
            nnod = nelem+1;
            vertex = zeros(4*nelem+1,2);  %%%%%% Ponerlo genérico
            
            for iNod = 1:nnod-1
                vertex(2*iNod-1,:)   = [0.2*z(iNod)/2 coord(iNod)];
                vertex(2*iNod,:) = [0.2*z(iNod)/2 coord(iNod+1)];
            end

            vertex(2*nelem+1:4*nelem,1) = - fliplr(vertex(1:2*nelem,1)')';
            vertex(2*nelem+1:4*nelem,2) = fliplr(vertex(1:2*nelem,2)')';
            vertex(end,:) = vertex(1,:);
            pgon = polyshape(vertex);
            tr = triangulation(pgon);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            geometryFromMesh(model,tnodes,telements);
            generateMesh(model,'GeometricOrder','linear','Hmax',0.1)
            coord = model.Mesh.Nodes';
            connec = model.Mesh.Elements';
            s.coord = coord;
            s.connec = connec;
            m = Mesh(s);

            figure(1);
            clf
            m.plot();
            grid on
            grid minor
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            

            xBar = obj.mesh.computeBaricenter();
            figure(2)
            subplot(2,2,[1 3]);plot(xBar,z)
            grid on
            grid minor
            title('Clamped-Clamped Column Profile','Interpreter', 'latex','FontSize',20, 'fontweight','b');
            xlabel('x','Interpreter', 'latex','fontsize',14,'fontweight','b');
            ylabel('A(x)','Interpreter', 'latex','fontsize',14,'fontweight','b');
        end

        function plotBucklingModes(obj,m1,m2)
            mod1 = m1;           
            mod2 = m2;
            coord = obj.mesh.coord;
            subplot(2,2,2); plot(coord,-mod1);
            grid on
            grid minor
            title('First Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
            subplot(2,2,4); plot(coord,-mod2);
            grid on
            grid minor
            title('Second Buckling Mode','Interpreter', 'latex','FontSize',14, 'fontweight','b')
        end

        function plotEigenvaluesAndIterations(obj)
            figure(3)
            clf
            hold on
            plot(1:obj.iter,obj.E1);
            plot(1:obj.iter,obj.E2);
            hold off
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Eigenvalues','Interpreter', 'latex','fontsize',18,'fontweight','b');
            axis([0 60 0 100]);
        end        
        
    end
    
end