Examples and test suite        
=======================

A number of small benchmark optimization test cases     
that illustrate the features available in the Null Space Optimizer package     
are available in the folder 
`examples <https://gitlab.com/florian.feppon/null-space-optimizer/-/tree/public-master/examples>`_.   
They can be executed as python modules from the command line with the `-m` option, e.g.:
    

.. code:: bash

   python -m nullspace_optimizer.examples.basic_examples.ex00_basic_problem 
    
or by directly executing the Python file in the example folder:
    
.. code:: bash  
    
   python ex00_basic_problem.py
    
These examples form a test suite that can be run from the ``nullspace_optimizer`` root folder:  

.. code:: bash

   pytest -vv

    
Basic examples  
--------------

.. subfigure:: ABC
   :layout-sm: A|B|C
   :gap: 3px
   :name: myfigure
   :class-grid: outline

   .. image:: img/ex02.png
      :alt: Image A

   .. image:: img/ex05.png
      :alt: Image B

   .. image:: img/ex06.png
      :alt: Image C

   Optimization trajectories obtained with the Null Space Optimizer for several     
   of the basic examples test cases.
    
The `basic_examples <https://gitlab.com/florian.feppon/null-space-optimizer/-/tree/public-master/examples/basic_examples>`_ folder contains minimal optimization test cases showing
how to implement     optimization problems and solving them with the Null Space
Optimizer.   The examples feature only two variables to test and illustrate
the method  with graphical visualizations of the optimization trajectories. 

.. dir:: examples/basic_examples


Advanced examples   
-----------------

The `advanced_optimizable_examples <https://gitlab.com/florian.feppon/null-space-optimizer/-/tree/public-master/examples/advanced_optimizable_examples>`_ folder gathers test cases using more advanced
features such as for automating    the    inclusion of bound constraints and
filtering the optimization variable.  For instance, the examples ``ex00_bound_constraints.py``,   
``ex01_bound_constraints_sparse.py`` and ``ex01_bound_constraints_sparse_cvxopt.py`` solve the following    
optimization problem:
            
.. math::     

 \begin{aligned}
     \min_{(x_{ij})_{1\leqslant i,j\leqslant n}}  & \quad J((x_{ij})_{1\leqslant i,j\leqslant n}):=-\frac{1}{2}
     \left[\sum_{i=1}^{n-1}\sum_{j=1}^{n}(x_{i+1,j}-x_{i,j})^{2}+\sum_{i=1}^{n}\sum_{j=1}^{n-1}(x_{i,j+1}-x_{i,j})^{2}\right]   
     \\ 
 s.t. & \quad  0\leqslant x_{ij}\leqslant 1 \text{ for all }1\leqslant i\leqslant n.
 \end{aligned}


.. figure:: img/checkerboard.png    
   :height: 250px
   :align: center
   :alt: Checkerboard benchmark problem result
    
   An optimized solution found with the Null Space Optimizer
    
   

    
.. dir::  examples/advanced_optimizable_examples

    
Topology optimization examples
------------------------------
    
.. subfigure:: ABC
   :layout-sm: A|B|C
   :gap: 3px
   :class-grid: outline

   .. image:: img/MBB.png
      :alt: Image A

   .. image:: img/multiple.png
      :alt: Image B
      :height: 100px

   .. image:: img/heat.png
      :alt: Image C
      :height: 100px

   Optimized designs computed with the Null Space Optimizer for several topology optimization test cases.

The `topopt_examples <https://gitlab.com/florian.feppon/null-space-optimizer/-/tree/public-master/examples/topopt_examples>`_   
folder contains examples code solving classical     
topology optimization test case with the density method.    
The folder includes classical examples in structural mechanics     
and heat conduction, solving single and multiple load cases and on structured and unstructured meshes with  
the Null Space Optimizer but also with MMA, IPOPT and the Optimality Criteria method.
The user is referred to the following reference  for an in-depth tutorial about the implementation  
of these examples:
    
.. pull-quote::

   Feppon F. *Density based topology optimization with the Null Space Optimizer: a
   tutorial and a comparison* (2023).   
   Submitted. HAL preprint `hal-04155507 <https://hal.archives-ouvertes.fr/hal-04155507/document>`_. 




    
.. dir::  examples/topopt_examples

