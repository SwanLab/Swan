function F=ellipse_system(AOA,del_x,del_y,del_ab)

    F=[cos(atan(del_ab(2)/del_ab(1))+AOA)*sqrt((del_ab(1))^2+(del_ab(2))^2)-del_x;
      sin(atan(del_ab(2)/del_ab(1))+AOA)*sqrt((del_ab(1)^2)+(del_ab(2)^2))-del_y];



end