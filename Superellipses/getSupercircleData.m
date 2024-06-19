function lhs = getSupercircleData(a,b,p,X,Y)
    lhs = abs((X) ./ a).^p + abs((Y) ./ b).^p;
end