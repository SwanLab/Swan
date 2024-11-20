classdef Data < handle

    properties (Access = public)
        nFeatures
        nSamples
        nLabels
        polyGrade
        Xtrain
        Ytrain       
        Xtest
        Ytest
        Ntest
    end

    properties (Access = private)
        X
        Y
        data
        fileName
        testRatio
        features
    end

    methods (Access = public)

        function obj = Data(cParams)            
            obj.init(cParams)
            obj.loadData();
            obj.buildModel();
            obj.splitdata();
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
                   obj.polyGrade = h.value;
                   obj.buildModel(obj.X,obj.polyGrade);
           end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName        = cParams.fileName;
            obj.testRatio       = cParams.testRatio;
            obj.polyGrade       = cParams.polynomialOrder;
            obj.features        = cParams.features;
        end

        function loadData(obj)
            f = obj.fileName;
            obj.data = readmatrix(f);
       %     fprintf('Features to be used (1:%d):',(size(obj.data,2)-1))
       %     feat = input(' ');
            xfeat = obj.features;
            x = obj.data(:, xfeat);
            yfeat = xfeat(end)+1:size(obj.data,2);
            y = obj.data(:, yfeat);
            obj.X = x;
            obj.Y = y;
            obj.nLabels   = size(obj.Y,2);                        
            obj.nFeatures = size(obj.X,2);
            obj.nSamples  = size(obj.X,1);

            % IDENTIFIER
            % ydata = obj.data(:, end);
            % y = zeros(length(ydata),max(ydata));

%             ydata = obj.data(:, feat);
%             y = zeros(length(ydata),width(ydata));
%             
%             u = unique(ydata);
%             for i=1:length(ydata)
%                 for j = 1:length(u)
%                     if ydata(i) == u(j)
%                         y(i,j) = 1;
%                     end
%                 end
%             end
%             
%             obj.X = (x-min(x,[],1))./(max(x,[],1)-min(x,[],1)+10^(-10));
%             % obj.Y = y;
%             obj.Y = obj.X;
        end

        function Xful = buildModel(obj)
            % Generation of all the possible exponent combinations
            exponents = [];
            for tD = obj.polyGrade:-1:1
                newDeg    = obj.generateExponents(tD);
                exponents = [exponents; newDeg];
            end
            exponents = flip(exponents,1);
            
            % Initialization of the output matrix
            Xful = zeros(obj.nSamples, size(exponents,1));
            
            % Double loop to cover all the possible polynomials
            for i = 1:size(exponents,1)
                auxTerm = ones(obj.nSamples,1);
                for j = 1:obj.nFeatures
                    auxTerm = auxTerm.*(obj.X(:,j).^exponents(i,j));
                end
                Xful(:,i) = auxTerm;
            end
            obj.X = Xful;
        end
        
        function exponents = generateExponents(obj,targetDeg)
            % Initialization of parameters
            exponents = [];
            currentExponents = zeros(1, obj.nFeatures);
            initialIndex = 1;
            
            % Calculation of the possible exponents for the target degree
            exponents = obj.generateExponentsRecursive(targetDeg,initialIndex,currentExponents,exponents);
        end
        
        function exponents = generateExponentsRecursive(obj,targetDeg,currentIndex,currentExponents,exponents)
            % Assignation of exponents for the base case
            if currentIndex == obj.nFeatures
                currentExponents(currentIndex) = targetDeg;
                exponents = [exponents; currentExponents];
            else
                % Recursion to search for the possibile combinations which
                % sum the polynomial degree target
                for i = 0:targetDeg
                    currentExponents(currentIndex) = i;
                    exponents = obj.generateExponentsRecursive(targetDeg - i, currentIndex + 1, currentExponents, exponents);
                end
            end
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