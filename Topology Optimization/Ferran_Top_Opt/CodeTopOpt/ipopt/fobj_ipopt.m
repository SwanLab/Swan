function [ f ] = fobj_ipopt( x,fobj )

global post_info

 if ~isequal(post_info.xold_fobj,x)
     [post_info.fobj,post_info.gradient_fobj] = fobj(x);
     post_info.xold_fobj = x;
 end
 f = post_info.fobj;

end

