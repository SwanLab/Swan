function lhs = getSupercircleData(a,b,p,X,Y)
    lhs = abs((X-0.5) ./ a).^p + abs((Y-0.5) ./ b).^p;
end