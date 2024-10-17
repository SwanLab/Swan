clear all
close all

correct_margin = false;
k=1;
h=1;
margin = k*(10^(-h));

while correct_margin == false


        if k<1
            h=h+1;
            k=9;
            disp('HOLAAAAAA')
        else
            k=k-0.5;
        end
        margin = k*(10^(-h));

    
    disp(margin)

end
if size(bMesh.coord,1)==0 || size(bMesh.connec,1)==0
    disp('MARGE NO TROBAT')
    return
end







