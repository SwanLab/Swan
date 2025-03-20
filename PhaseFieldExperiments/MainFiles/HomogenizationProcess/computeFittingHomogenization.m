function [fun,dfun,ddfun] = computeFittingHomogenization(cParams,nonZeroComp,type)
    x = cParams.phi;
    x = reshape(x,length(x),[]);
    y = cParams.mat;
    
    fun   = cell(3,3);
    dfun  = cell(3,3);
    ddfun = cell(3,3);
    for i=1:3
        for j=1:3
            if i==1 && j==1 && nonZeroComp
                switch type
                    case 1
                        coeffs = polyfit(x,squeeze(y(i,j,:)),9);
                    case 2
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                    case 3
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[squeeze(y(i,j,end))]);
                    case 4
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),squeeze(y(i,j,end))]);
                end
            else
                switch type
                    case 1
                        coeffs = polyfit(x,squeeze(y(i,j,:)),9);
                    case 2
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                    case 3
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[0]);
                    case 4
                        coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),0]);
                end
            end
            fun{i,j} = poly2sym(coeffs);
            dfun{i,j} = diff(fun{i,j});
            ddfun{i,j} = diff(dfun{i,j});
        end
    end
 end