clear;

acc = "Forward";

switch acc

    case "Forward"

        x = ValDerForward(2,[1 0]);
        y = ValDerForward(0.5,[0 1]);
        f = 0.5 * u * A * u;

        vec = f.double;

        disp("The calculation of the function is " + vec(1));
        disp("The calculation of the derivative respect x of the function is " + vec(2));
        disp("The calculation of the derivative respect y of the function is " + vec(3));

    case "Reverse"

        x = ValDerReverse(2);
        y = ValDerReverse(0.5);
        f = cos(x)+x*y^2+x^2;

    otherwise

        disp("Method is not valid.");
        return

end


