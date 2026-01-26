function [out] = sym2str(sy)
%updated:  02-03-2009
%author: Marty Lawson
%
%converts symbolic variables to a matlab equation string insuring that
%only array opps are used.  
%Symbolic arrays are converted to linear Cell arrays of strings. 
%This function is especially usefull when combined with the eval() command.  
%Also, converts maple atan function to matlab atan2 and converts 
%maple "array([[a,b],[c,d]])" notation to matlab "[a,b;c,d]" notation.  
%
%Note: eval() of a matrix of functions only works if all the input 
%variables have single values.  i.e. vectors and arrays won't work.
%
%Note2: eval() does not work on Cell arrays directly.  Use "Cell_array{index}"
%inside of the eval() to keep eval() happy
%
%	EXAMPLE:
%
% X_t = dsolve('5*D2x+6*Dx+3*x = 10*sin(10*t)','x(0)=0','Dx(0)=3'); %solution is a symbolic function, X(t)
% X_t_str = sym2str(X_t);							                %convert from symbolic type to char type using array operations
% t = 0:.01:20;									                    %make "t" an array containing the time range of interest for X(t)
% X_t_vec = eval(X_t_str);							                %see "help eval" for details.
% plot(t,X_t_vec)									                %plot the results
% grid on										                    %make the plot look nice
% xlabel('time [radians]')
% ylabel('amplitude [-]')
% 
    sy = sym(sy); %insure input is symbolic
    siz = prod(size(sy));   %find the number of elements in "sy"
    for i = 1:siz   %dump it into a cell array with the same number of elements
        in{i} = char(sy(i)); %convert to char
        in{i} = strrep(in{i},'^','.^');%insure that all martix opps are array opps
        in{i} = strrep(in{i},'*','.*');
        in{i} = strrep(in{i},'/','./');
        in{i} = strrep(in{i},'atan','atan2'); %fix the atan function
        in{i} = strrep(in{i},'array([[','['); %clean up any maple array notation
        in{i} = strrep(in{i},'],[',';');
        in{i} = strrep(in{i},']])',']');
    end
    if siz == 1
        in = char(in);      %revert back to a 'char' array for single answers
    end
    out = in;