function Xful = buildModel(x,d)
    x1 = x(:,1);
    x2 = x(:,2);
    cont = 1;
    for g = 1:d
        for a = 0:g
            Xful(:,cont) = x2.^(a).*x1.^(g-a);
            cont = cont+1;
        end
    end
end