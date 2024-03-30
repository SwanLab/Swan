clear;

acc = "Forward";

switch acc

    case "Forward"

        x = ValGradForward(2,[1 0]);
        y = ValGradForward(0.5,[0 1]);
        f = cos(x)+x*y^2+x^2;

        vec = f.double;

        disp("The calculation of the function is " + vec(1));
        disp("The calculation of the Gradivative respect x of the function is " + vec(2));
        disp("The calculation of the Gradivative respect y of the function is " + vec(3));

    case "Reverse"

        x = ValGradReverse(2);
        y = ValGradReverse(0.5);
        f = cos(x)+x*y^2+x^2;

    otherwise

        disp("Method is not valid.");
        return

end


