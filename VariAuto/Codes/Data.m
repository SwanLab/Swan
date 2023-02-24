classdef Data < handle

    properties (Access = public)
        nData
        nFeatures
        nLabels
        data
        Xtrain
        Ytrain       
        Xtest
        Ytest
        polyGrade
    end

    properties (Access = private)
        testRatio
        X
        Y
    end

    methods (Access = public)
        function self = Data(File_Name,TP,d)
            self.init(File_Name,TP,d)
        end

        function plotdata(self,i,j)
            gscatter(self.Xtrain(:,i),self.Xtrain(:,j),self.Ytrain,'rgbcmyk','*')
            xlabel(['X',num2str(i)]);
            ylabel(['X',num2str(j)]);
        end

        function plotCorrRow(self,k)
            x = self.data(:,1:end-1);
            nf = size(x,2);
            for i = 1:nf
                nexttile((i-1)*nf+k)
                if i == k
                    histogram(x(:,i))
                else
                    gscatter(self.data(:,k),self.data(:,i),self.Y,'rgbcmyk','*')
                end
            end
        end

        function plotCorrMatrix(self)
            x = self.data(:,1:end-1);
            nf = size(x,2);
            figure            
            t = tiledlayout(nf,nf,'TileSpacing','Compact'); 
            title(t,'Features correlation matrix');
            for i = 1:nf
                self.plotCorrRow(i);
            end
        end

        function updateHyperparameter(self,h)
           switch h.type
               case 'testRatio'
                   self.testRatio = h.value;
                   self.splitdata()
               case 'polyGrade'
                   self.polyGrade = h.value;
                   self.computefullvars(self.X,self.polyGrade);
           end
        end
    end

    methods (Access = private)

        function init(self,File_Name,TP,d)
            self.testRatio = TP;
            self.polyGrade = d;
            self.loadData(File_Name);
            self.computefullvars(self.X,d)
            self.splitdata()
            self.nData = size(self.data,1);
            self.nLabels = size(self.Ytrain,2);                        
            self.nFeatures = size(self.Xtrain,2);
        end

        function loadData(self,FN)
            f = fullfile('../Datasets/', FN);
            self.data = load(f);
            fprintf('Features to be used (1:%d):',(size(self.data,2)-1))
            feat = input(' ');
            x = self.data(:, feat);
            ydata = self.data(:, end);
            y = zeros(length(ydata),max(ydata));
            u = unique(ydata);
            for i=1:length(ydata)
                for j = 1:length(u)
                    if ydata(i) == u(j)
                        y(i,j) = 1;
                    end
                end
            end
            self.X = (x-min(x,[],1))./(max(x,[],1)-min(x,[],1)+10^(-10));
            self.Y = y;
%             self.X = x; 
%             self.Y = ydata;
        end
        
        function computefullvars(self,x,d)
            self.X = self.X;
%             self.X = buildModel(x,d);
            self.Y = self.Y;
        end 

        function splitdata(self)
            nD = size(self.data,1);
            TP = self.testRatio;
            r = randperm(nD);
            ntest = round(TP/100*nD);
            ntrain = nD - ntest;
            self.Xtrain = self.X(r(1:ntrain),:);
            self.Xtest = self.X(r((ntrain + 1):end),:);
            self.Ytrain = self.Y(r(1:ntrain),:);
            self.Ytest = self.Y(r((ntrain + 1):end),:);
        end
    end
end