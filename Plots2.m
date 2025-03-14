matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');
matType{6} = load('HorizontalCrackMicroDamageArea.mat');
matType{6}.mat = matType{6}.mat*210;
matType{6}.phi = matType{6}.holeParam{1};

 for i=1:length(matType)
     [funMat(:,:,i),dfunMat(:,:,i),ddfunMat(:,:,i)] = computeFunctionsAndDerivatives(matType{i});
 end

 function [fun,dfun,ddfun] = computeFunctionsAndDerivatives(cParams)
    x = cParams.phi;
    x = reshape(x,length(x),[]);
    y = cParams.mat;
    
    fun   = cell(3,3);
    dfun  = cell(3,3);
    ddfun = cell(3,3);
    for i=1:3
        for j=1:3
            coeffs = polyfit(x,squeeze(y(i,j,:)),9);
            %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
            %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),0]);
            %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[0]);
            fun{i,j} = poly2sym(coeffs);
            dfun{i,j} = diff(fun{i,j});
            ddfun{i,j} = diff(dfun{i,j});
        end
    end
end