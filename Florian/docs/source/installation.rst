Installation
============

You can install the Null Space Optimizer    
with pip or by cloning the sources from the gitlab repository. 
    

Installation Requirements   
-------------------------
    
For installing Null Space Optimizer, a python executable (version >=3.6)    
needs to be installed and available from the command line.    

If you don't install Null Space Optimizer with ``pip``, you will also need to manually install the 
following python packages:
    
- `numpy <https://numpy.org/>`_ (>=1.23.4)
- `matplotlib <https://matplotlib.org/>`_ (>=3.6.2)
- `scipy <https://scipy.org/>`_ (>=1.9.3) 
- `cvxopt <https://cvxopt.org/install/index.html/>`_ (>=1.3.0)    
- `osqp <https://osqp.org/docs/get_started/python.html>`_ (>=0.6.2)
- `sympy <https://docs.sympy.org/latest/>`_ (>=1.3.0)
- `qpalm <https://github.com/kul-optec/QPALM>`_ (>=1.2.2)
- `colored <https://dslackw.gitlab.io/colored/>`_ (>=1.4.4)


Installation with pip
---------------------
    
This is the easiest way to get the latest stable version and to install python dependencies.

.. code:: bash

   pip install nullspace_optimizer
    
    
Installation from the gitlab repository
---------------------------------------

Clone the repository using the command line prompt:

.. code:: console

   git clone https://gitlab.com/florian.feppon/null-space-optimizer.git

Then add the ``nullspace_optimizer`` folder to your ``$PYTHONPATH``:  
    
.. code-block:: bash  
   :caption: .bashrc    
    
   export PYTHONPATH="$PYTHONPATH:/path/to/nullspace_optimizer"

or install the module locally with ``pip``:

.. code:: bash

   pip install -e /path/to/nullspace_optimizer

where ``/path/to/nullspace_optimizer`` is the directory where Null Space Optimizer has
been cloned.
    

Installation requirements for running topology optimization examples    
--------------------------------------------------------------------
    
If you want to run the topology optimization examples from the ``examples/topopt_examples``, you    
also need to install the following dependencies:    
    
- `IPOPT <https://coin-or.github.io/Ipopt/INSTALL.html>`_   
  and its Python wrapper `cyipopt <https://cyipopt.readthedocs.io/en/stable/install.html#from-source>`_. Building from  
  source is recommended.
- `pymedit <https://gitlab.com/florian.feppon/pymedit>`_
- `FreeFEM <https://doc.freefem.org/introduction/installation.html>`_ and   
  its Python interface `PyFreeFEM <https://pyfreefem.readthedocs.io/en/latest/installation.html>`_.     
  Building FreeFEM from source is also recommended. An installation 
  recipe for Linux systems is   
  described `here <https://people.cs.kuleuven.be/~florian.feppon/topopt_course/install_freefem.html>`_.
    
