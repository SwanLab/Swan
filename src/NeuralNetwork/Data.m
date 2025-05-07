classdef Data < handle

    properties (Access = public)
        nFeatures
        nSamples
        nLabels        
        Xtrain
        Ytrain       
        Xtest
        Ytest
        Ntest
    end

    properties (Access = private)
        X
        Y
        polynomialOrder
        data
        fileName
        testRatio
        xFeatures
        yFeatures
        k
    end

    methods (Access = public)

        function obj = Data(cParams)            
            obj.init(cParams)
            obj.loadData();
            obj.buildModel();
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
            obj.xFeatures       = cParams.xFeatures;
            obj.yFeatures       = cParams.yFeatures;
            obj.k               = cParams.k;
        end

        function loadData(obj)

            obj.data = load(obj.fileName).data;

            % Change: incorporate features to use in cParams vs propmpting
            % user though terminal
            x = obj.data(:, obj.xFeatures);
            y = obj.data(:, obj.yFeatures);

            obj.X = x;
            obj.Y = y;

        end
        

        function Xful = buildModel(obj)
            x    = obj.X;
            Xful = x;
            d    = obj.polynomialOrder;

            if d > 1
                
                for i = 2:d
                    tempX = x.^i;
                    Xful = cat(2,Xful,tempX);
                end
            
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
            ntest = round(TP/100*nD);
            TestStart = 1 + ntest * (obj.k - 1.0);
            TestEnd   = min(ntest * obj.k,nD); 
            
            testIdx = TestStart:TestEnd;
            allIdx  = 1:nD;
            trainIdx = setdiff(allIdx, testIdx);
            obj.Xtrain = obj.X(trainIdx, :);
            obj.Ytrain = obj.Y(trainIdx, :);
            obj.Xtest  = obj.X(testIdx, :);
            obj.Ytest  = obj.Y(testIdx, :);
            obj.Ntest = ntest;
        end
    end

    methods (Access = public, Static)

        function [NewData] = preProcessDataset(data)

            nD = size(data,1);
            r  = randperm(nD);
            NewData = data(r(1:nD),:);

            maxValue    = max(NewData(:,end));
            minValue    = min(NewData(:,end));
            NewData(:,end)  = (NewData(:,end) - minValue) / (maxValue - minValue);

        end

    end
end