function data = mginit(pv, hmax, nref)
    [p, t, e] = pmesh(pv, hmax, 1);
    data(1).p = p; 
    data(1).t = t; 
    data(1).e = e;
    
    for i = 1:nref+1
        [data(i).T] = interpolation(data(i).p,data(i).t);
        data(i+1).p = data(i).T*data(i).p; 
        data(i+1).t = delaunayn(data(i+1).p);
        data(i+1).e = boundary_nodes(data(i+1).t);
    end   
    
    for i = 1:nref+1
        data(i).R = restriction(data(i).T);          
    end
    
    [~,A,b] = fempoi(data(nref+1).p,data(nref+1).t,data(nref+1).e);
    data(nref + 1).A = A;
    data(nref + 1).b = b; 
    
    for j = nref:-1:1
        data(j).A = data(j).R * data(j+1).A * data(j).T;   
    end
    
end
    
    
    