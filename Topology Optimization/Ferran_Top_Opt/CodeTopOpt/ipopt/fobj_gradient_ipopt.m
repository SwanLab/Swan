function [ g ] = fobj_gradient_ipopt( x,fobj )

global post_info

 if ~isequal(post_info.xold_fobj,x)
     [post_info.fobj,post_info.gradient_fobj] = fobj(x);
     post_info.xold_fobj = x;
 end
 g = post_info.gradient_fobj;

end

