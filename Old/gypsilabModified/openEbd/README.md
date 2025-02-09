# EBD : Efficient Bessel Decomposition toolbox (matlab)
#### Fast convolution with radial kernels in 2D


This toolbox implements the algorithm described in "Discrete convolution in $\mathbb{R}^2$  with radial kernels using non-uniform fast Fourier transform with non-equispaced frequencies", written by Martin Averseng, and submitted to the journal Numerical Algorithms in 2018. 

To test it, you can directly run Demo.m, or DemoGrad. Here follows a more detailed description of the algorithm and a tutorial. The method is designed to compute fast approximations of vectors $$q$$ which entries are given by 
<a href="https://www.codecogs.com/eqnedit.php?latex=$$q_k&space;=&space;\sum_{l=1}^{N_y}&space;G(X_k&space;-&space;Y_l)&space;f_l,&space;\quad&space;k&space;=&space;1,&space;\cdots,&space;N_x$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$q_k&space;=&space;\sum_{l=1}^{N_y}&space;G(X_k&space;-&space;Y_l)&space;f_l,&space;\quad&space;k&space;=&space;1,&space;\cdots,&space;N_x$$" title="$$q_k = \sum_{l=1}^{N_y} G(X_k - Y_l) f_l, \quad k = 1, \cdots, N_x$$" /></a>

or

<a href="https://www.codecogs.com/eqnedit.php?latex=$$q_k&space;=&space;\sum_{l=1}^{N_y}&space;\nabla&space;G(X_k&space;-&space;Y_l)&space;f_l,&space;\quad&space;k&space;=&space;1,&space;\cdots,&space;N_x$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$q_k&space;=&space;\sum_{l=1}^{N_y}&space;\nabla&space;G(X_k&space;-&space;Y_l)&space;f_l,&space;\quad&space;k&space;=&space;1,&space;\cdots,&space;N_x$$" title="$$q_k = \sum_{l=1}^{N_y} \nabla G(X_k - Y_l) f_l, \quad k = 1, \cdots, N_x$$" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=$$G$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$G$" title="$G$" /></a> is a radial function i.e. <a href="https://www.codecogs.com/eqnedit.php?latex=$$G(x)&space;=&space;g(\lvert&space;x\rvert)$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$G(x)&space;=&space;g(\lvert&space;x\rvert)$$" title="$$G(x) = g(\lvert x\rvert)$$" /></a> for some function <a href="https://www.codecogs.com/eqnedit.php?latex=$g$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$g$$" title="$$g$$" /></a>,<a href="https://www.codecogs.com/eqnedit.php?latex=$$X$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$X$$" title="$$X$$" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=$$Y$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$Y$$" title="$$Y$$" /></a>  are two clouds of<a href="https://www.codecogs.com/eqnedit.php?latex=$$N_x$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$N_x$$" title="$$N_x$$" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=$$N_y$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$N_y$$" title="$$N_y$$" /></a> points in <a href="https://www.codecogs.com/eqnedit.php?latex=$$\mathbb{R}^2$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\mathbb{R}^2$$" title="$$\mathbb{R}^2$$" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=$$f$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$f$$" title="$$f$$" /></a> is a complex vector. 

### 1°) Description of the algorithm

The method first decomposes G in finite Bessel series 

<a href="https://www.codecogs.com/eqnedit.php?latex=$$G(r)&space;=&space;\sum_{p&space;=&space;1}^P&space;\alpha_p&space;J_0(\rho_p&space;r)&space;$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$G(r)&space;=&space;\sum_{p&space;=&space;1}^P&space;\alpha_p&space;J_0(\rho_p&space;r)&space;$$" title="$$G(r) = \sum_{p = 1}^P \alpha_p J_0(\rho_p r) $$" /></a>

where J0 is the Bessel function of first kind, rho is the sequence of its positive zeros, and alpha are called the EBD coefficients of G. The coefficients are chosen as the minimizers of the Sobolev H^10 norm error in this approximation on a ring rmin < r < rmax where rmin is a cutoff parameter and rmax is the greatest distance occurring between two points Xk and Yl. The method takes the paramter a = rmin/rmax as an input. 

Then, each J0(rhop r) is approximated by 

<a href="https://www.codecogs.com/eqnedit.php?latex=$$J0(\rho_p&space;|x|)&space;=&space;\frac{1}{M_p}\sum_{m&space;=&space;1}^{M_p}&space;e^{i&space;\rho_p&space;\xi_{m}^p&space;\cdot&space;x}$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$J0(\rho_p&space;|x|)&space;=&space;\frac{1}{M_p}\sum_{m&space;=&space;1}^{M_p}&space;e^{i&space;\rho_p&space;\xi_{m}^p&space;\cdot&space;x}$$" title="$$J0(\rho_p |x|) = \frac{1}{M_p}\sum_{m = 1}^{M_p} e^{i \rho_p \xi_{m}^p \cdot x}$$" /></a>

where xi^pm = exp(i 2mpi/Mp). This is the trapezoidal rule applied to the formula 

<a href="https://www.codecogs.com/eqnedit.php?latex=$$J_0(|x|)&space;=&space;\int_{\partial&space;B}&space;e^{i&space;x&space;\cdot&space;\xi}&space;d\xi$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$J_0(|x|)&space;=&space;\int_{\partial&space;B}&space;e^{i&space;x&space;\cdot&space;\xi}&space;d\xi$$" title="$$J_0(|x|) = \int_{\partial B} e^{i x \cdot \xi} d\xi$$" /></a>

where the integration takes place on the boundary of the unit disk B in R2.

Combining these two steps, we obtain an approximation for G of the form 

<a href="https://www.codecogs.com/eqnedit.php?latex=$$G(x)&space;\approx&space;\sum_{\nu&space;=&space;1}^{N_\xi}&space;\hat{\omega}_\nu&space;e^{i&space;\xi_\nu&space;\cdot&space;x}$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$G(x)&space;\approx&space;\sum_{\nu&space;=&space;1}^{N_\xi}&space;\hat{\omega}_\nu&space;e^{i&space;\xi_\nu&space;\cdot&space;x}$$" title="$$G(x) \approx \sum_{\nu = 1}^{N_\xi} \hat{\omega}_\nu e^{i \xi_\nu \cdot x}$$" /></a>

valid for |x| > rmin. If this is replaced in the expression of qk, we see that q can be approximated by non-uniform Fourier transform for any vector f. The interactions |Xk - Yl| < rmin where the approximation is not valid are corrected by a sparse matrix product. 


The code is in Matlab language. The implementation of the NUFFT is borrowed from Leslie Greengard, June-Yub Lee and Zydrunas Gimbutas (see license file in the libGgNufft2D folder). The ideas come from a similar method in 3D called Sparse Cardinal Sine Decomposition, developped by François Alouges and Matthieu Aussal, also published in Numerical Algorithms. 

### 2°) Tutorial

#### A. Simple EBD
- Create the arrays X and Y of sizes Nx x  2 and Ny x 2 (points in R^2). 
- Create a kernel by calling 
```
G = Kernel(fun,der);
```

where fun is an anonymous function of your choice and der is an anonymous function 
repesenting the derivative of fun. For example, 
```
G = Kernel(@(x)(1./x),@(x)(-1./x.^2));
```

Be aware that the anonymous functions provided in argument must accept arrays as input. For some classical kernels, the methods have been optimized. You can create one of those special kernels using the Kernel library:
```
G = LogKernel;  % represents G(x) = log(x);
G = ThinPlate(a,b) % represents G(x) = a*x^2*log(b*x);
```

(and others, see folder Kernels).

- Define the tolerance in the error of approximation. The method guarantees that 
the Bessel decomposition of G in the ring is accurate at the tolerance level. 
This implies that the maximal error on an entry of qk is `tol*norm(q,1)`.  

- Define the parameter a (the ratio between rmin and rmax.) It is roughly the proportion of interactions that will be computed exactly. There is an optimal value of a for which the evaluation of the convolution is the fastest, but it is not possible to know it in advance. However, when X and Y are uniformly distributed on a disk, the optimal a is of the order <a href="https://www.codecogs.com/eqnedit.php?latex=$$\frac{1}{(N_x&space;N_y)^{1/4}}$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\frac{1}{(N_x&space;N_y)^{1/4}}$$" title="$$\frac{1}{(N_x N_y)^{1/4}}$$" /></a>, and if they are uniformly distributed on a curve, the optimal a is of the order <a href="https://www.codecogs.com/eqnedit.php?latex=$$\frac{1}{(N_x&space;N_y)^{1/3}}$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\frac{1}{(N_x&space;N_y)^{1/3}}$$" title="$$\frac{1}{(N_x N_y)^{1/3}}$$" /></a>. If a is large, a lot of interactions are computed exactly while the Bessel decomposition will have only a few terms. If a is small, the opposite will happen. 

- You can now call
```
[onlineEBD, rq, loc] = offlineEBD(G,X,Y,a,tol);
```

The variable `onlineEBD` contains a handle function. If f is a vector with length equal to size(Y,1)
```
q = onlineEBD(f);
```
returns the approximation of the convolution. The variable `rq` contains a RadialQuadrature object. You can visualize a representation of the radial approximation with
```
 rq.show();
```
and also check all its properties (including coefficients, frenquencies used, accuracy, ...). Finally, `loc` is the local correction matrix (used inside onlineEBD). 

### B. Derivative EBD 

If you rather want to compute the vectors 

<a href="https://www.codecogs.com/eqnedit.php?latex=$$[q_1(k),q_2(k)]&space;=&space;\sum_{l=1}^{N_y}&space;\nabla{G}(Y_l&space;-&space;X_k)f_l$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$[q_1(k),q_2(k)]&space;=&space;\sum_{l=1}^{N_y}&space;\nabla{G}(Y_l&space;-&space;X_k)f_l$$" title="$$[q_1(k),q_2(k)] = \sum_{l=1}^{N_y} \nabla{G}(Y_l - X_k)f_l$$" /></a>

where for a point x in <a href="https://www.codecogs.com/eqnedit.php?latex=$$\mathbb{R}^2$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\mathbb{R}^2$$" title="$$\mathbb{R}^2$$" /></a>,  

<a href="https://www.codecogs.com/eqnedit.php?latex=$$\nabla&space;G(x)&space;=&space;G'(|x|)\frac{&space;x}{|x|}$$," target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\nabla&space;G(x)&space;=&space;G'(|x|)\frac{&space;x}{|x|}$$," title="$$\nabla G(x) = G'(|x|)\frac{ x}{|x|}$$," /></a>

you may call `offline_dEBD`instead of `offlineEBD`. The syntax is

```
[MVx,MVy, rq, loc] = offlineEBD(G,X,Y,a,tol);
```
where, for example, `MVx`contains the handle function such that 
```
q1 = MVx(f)
```
is the component q1 of the convolution with \nabla G. 

    