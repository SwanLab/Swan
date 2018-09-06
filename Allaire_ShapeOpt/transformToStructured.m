function [mat] = transformToStructured(vect)
load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));

if length(b1) == length(vect)
    b = b1;
else 
    b = b0;
end

if length(dim) == 2
    for n = 1:length(vect)
        mat(b(n,1),b(n,2)) = vect(n);
    end
else
    for n = 1:length(vect)
        mat(b(n,1),b(n,2),b(n,3)) = vect(n);
    end 
end
end

