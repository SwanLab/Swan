classdef Data < handle

    properties (Access = public)
        nFeatures
        nLabels
        
        Xtrain
        Ytrain       
        Xtest
        Ytest
        Ntest
    end

    properties (Access = private)
        polynomialOrder
        X
        Y
        data
        fileName
        testRatio
    end

    methods (Access = public)

        function obj = Data(cParams)            
            obj.init(cParams)
            obj.loadData();
            obj.splitdata()
            obj.nLabels   = size(obj.Ytrain,2);                        
            obj.nFeatures = size(obj.Xtrain,2);            
        end

        function plotdata(self,i,j)
            gscatter(self.Xtrain(:,i),self.Xtrain(:,j),self.Ytrain,'rgbcmyk','*')
            xlabel(['X',num2str(i)]);
            ylabel(['X',num2str(j)]);
        end

        function plotCorrRow(obj,k)
            x = obj.data(:,1:end-1);
            nf = size(x,2);
            for i = 1:nf
                nexttile((i-1)*nf+k)
                if i == k
                    histogram(x(:,i))
                else
                    gscatter(obj.data(:,k),obj.data(:,i),obj.Y,'rgbcmyk','*')
                end
            end
        end

        function plotCorrMatrix(obj)
            x = obj.data(:,1:end-1);
            nf = size(x,2);
            figure            
            t = tiledlayout(nf,nf,'TileSpacing','Compact'); 
            title(t,'Features correlation matrix');
            for i = 1:nf
                obj.plotCorrRow(i);
            end
        end

        function updateHyperparameter(obj,h)
           switch h.type
               case 'testRatio'
                   obj.testRatio = h.value;
                   obj.splitdata()
               case 'polyGrade'
                   obj.polynomialOrder = h.value;
                   obj.buildModel(obj.X,obj.polynomialOrder);
           end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName        = cParams.fileName;
            obj.testRatio       = cParams.testRatio;
            obj.polynomialOrder = cParams.polynomialOrder;
        end

        function loadData(obj)
            f = fullfile('../Datasets/',obj.fileName);
            obj.data = load(f);
            fprintf('Features to be used (1:%d):',(size(obj.data,2)-1))
            feat = input(' ');
            x = obj.data(:, feat);

            % IDENTIFIER
            % ydata = obj.data(:, end);
            % y = zeros(length(ydata),max(ydata));

            ydata = obj.data(:, feat);
            y = zeros(length(ydata),width(ydata));
            
            u = unique(ydata);
            for i=1:length(ydata)
                for j = 1:length(u)
                    if ydata(i) == u(j)
                        y(i,j) = 1;
                    end
                end
            end
            
            obj.X = (x-min(x,[],1))./(max(x,[],1)-min(x,[],1)+10^(-10));
            % obj.Y = y;
            obj.Y = obj.X;
        end
        

        function Xful = buildModel(obj)
            x  = obj.X;
            d  = obj.polynomialOrder;
            x1 = x(:,1);
            x2 = x(:,2);
            cont = 1;
            for g = 1:d
                for a = 0:g
                    Xful(:,cont) = x2.^(a).*x1.^(g-a);
                    cont = cont+1;
                end
            end
            obj.X = Xful;
        end

        function splitdata(obj)
            nD = size(obj.data,1);
            TP = obj.testRatio;
            r = randperm(nD);
            ntest = round(TP/100*nD);
            ntrain = nD - ntest;
            obj.Xtrain = obj.X(r(1:ntrain),:);
            obj.Xtest  = obj.X(r((ntrain + 1):end),:);
            obj.Ytrain = obj.Y(r(1:ntrain),:);
            obj.Ytest  = obj.Y(r((ntrain + 1):end),:);
            obj.Ntest = ntest;
        end
    end
end