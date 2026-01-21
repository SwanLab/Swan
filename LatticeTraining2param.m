classdef LatticeTraining2param < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        h1
        h2
        tFrame
        Nr
        Ntheta
        nFrame
        connec
        CoarseOrder
        Length
    end

    methods (Access = public)

        function obj = LatticeTraining2param()
            close all
            obj.init()
            for i = 1:length(obj.h1)
                for j=1:length(obj.h2)
                    m = obj.createReferenceMesh(obj.h1(i),obj.h2(j));
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
                    EIFEoper.h1 = obj.h1(i);
                    EIFEoper.h2 = obj.h2(j);
                    %                 filePath = ['/home/raul/Documents/GitHub/EPFL/test/data_' num2str(obj.r(i), '%.3f') '.mat'];
                    filePath = ['./EPFL/dataLattice2param/data_h1_' num2str(obj.h1(i), '%.3f') '_h2_' num2str(obj.h2(j), '%.3f') '.mat'];
                    save(filePath,'EIFEoper')
                end
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
            obj.h1    = 0.01:0.1:1.95;
            obj.h2    = 0.01:0.1:1.95;
            obj.tFrame = 0.02;
            obj.Nr=10;
            obj.Ntheta=10;
            obj.nFrame = 1;
            obj.CoarseOrder = 1;
            obj.Length = 2; 
        end

        function printdisplacements(obj,Usbd,mesh,ind)
            for i = 1:size(Usbd,2)
                v = Usbd(:,i);
                EIFEMtesting.plotSolution(v,mesh,1,ind,i,[])
            end
        end

        function mesh = createReferenceMesh(obj,h1,h2)
            % X mesh
%             mesh =  meshLattice2param(obj.Length/2,obj.tFrame,h1,h2,obj.Nr,obj.Ntheta);
            % + mesh
            mesh = meshCrossLattice(obj.Length/2,obj.tFrame,h1,h2,obj.Nr,obj.Ntheta,obj.nFrame);

        end




    end

end