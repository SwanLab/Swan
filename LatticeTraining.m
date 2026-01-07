classdef LatticeTraining < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        h
        xmin
        xmax
        ymin
        ymax
        cy
        cx
        Nr
        Ntheta
        connec
        CoarseOrder
        Length
    end

    methods (Access = public)

        function obj = LatticeTraining()
            close all
            obj.init()

            for i=1:length(obj.h)
                m = obj.createReferenceMesh(obj.h(i));
%                  m = obj.meshFromPython(obj.r(i));
                if i==1
                    obj.connec = m.connec;
                end

                if norm(m.connec-obj.connec)>0
                    break
                end

                %                 figure
                %                 m.plot()
                %                 m.plotAllNodes();
                %                 m = obj.createReferenceMesh();
                data = Training(m,obj.CoarseOrder);
%                                 obj.printdisplacements(data.uSbd,m,i)
                p = OfflineDataProcessor(data);
                EIFEoper = p.computeROMbasis();
                EIFEoper.U = EIFEoper.Udef + EIFEoper.Urb;
                EIFEoper.U = EIFEoper.U(:);
%                 EIFEoper.Kfine = data.LHSsbd;
                EIFEoper.snapshots = data.uSbd;
%                 filePath = ['/home/raul/Documents/GitHub/EPFL/test/data_' num2str(obj.r(i), '%.3f') '.mat'];
                filePath = ['./EPFL/dataLattice/data_' num2str(obj.h(i), '%.3f') '.mat'];
                save(filePath,'EIFEoper')
            end
        end

    end

    methods (Access = private)

        function init(obj)
%             N = 80;
%             % Interval bounds
%             a = 0.99;
%             b = 0.8;
%             % Index vector
%             i = 0:N;
%             % Cosine spacing formula
%             obj.r = (a + b)/2 + (b - a)/2 * cos(pi * (1 - i / N));
            obj.h    = 0.8:0.01:0.96;
            obj.xmin = -1;
            obj.xmax = 1;
            obj.ymin = -1;
            obj.ymax = 1;
            obj.cx = 0;
            obj.cy = 0;
            obj.Nr=10;
            obj.Ntheta=10;
            obj.CoarseOrder = 1;
            obj.Length = 1;
        end

        function printdisplacements(obj,Usbd,mesh,ind)
            for i = 1:size(Usbd,2)
                v = Usbd(:,i);
                EIFEMtesting.plotSolution(v,mesh,1,ind,i,[])
            end
        end

        function mesh = createReferenceMesh(obj,h)

            %UnitMesh better
            %             x1      = linspace(-1,1,50);
            %             x2      = linspace(-1,1,50);
            %             [xv,yv] = meshgrid(x1,x2);
            %             [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            %             s.coord  = V(:,1:2);
            %              s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)-[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)-[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,0];
            %             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
            %                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,0];
            %             s.connec = F;
            %             mesh = Mesh.create(s);
            m =  mesh_square_X_solid(obj.Length,h,obj.Nr,obj.Ntheta);
            delta= 1e-9;
            s.coord = m.coord;
%             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
%                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
%             s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
%                 s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,+delta];
%             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
%                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[+delta,-delta];
%             s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
%                 s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[+delta,+delta];
            s.connec=m.connec;
            mesh = Mesh.create(s);

        end

        

       
    end

end