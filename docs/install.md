
# Install & Use

  [Home](index.md) |
  [Install&Use](install.md) |
  [About](about.md)

## Installation Prerequisities

Our framework is written in [**Python**](https://www.python.org/) programming language. By default, **Python** is available on any Linux platform. The optimal assignment of **Inferred Scaffolding** to **Reference Scafolding** is performed by using [**IBM CPLEX**](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) optimizer.

We use the following scripts/libraries:

- **Fastaq** script suite from [(Hunt et al., *Genome Biology*, 2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r42). For convenience, the source of **Fastaq** is included in the repository
- [Biopython](http://biopython.org/) library
- Python API of CPLEX

### Install Biopython

To install Biopython, type the following commands in your prompt:

~~~bash
pip install numpy
pip install biopython
~~~

### Install Python API of CPLEX

**IBM CPLEX** optimizer is a very powerful tool used for solving a wide range of optimization problems (Linear Programming, etc). A free CPLEX academic license can be obtained [here](https://www.ibm.com/developerworks/community/blogs/jfp/entry/cplex_studio_in_ibm_academic_initiative). As per **CPLEX** documentation:

To install the **CPLEX-Python** modules on your system, use the script **setuy.py** located in **yourCplexhome/python/PLATFORM**. If you want to install the **CPLEX-Python** modules in a nondefault location, use the option **--home** to identify the installation directory. For example, to install the **CPLEX-Python** modules in the default location, use the following command from the command line:

~~~bash
python setup.py install
~~~

To install in the directory **yourPythonPackageshome/cplex**, use the following command from the command line:

~~~bash
python setup.py install --home yourPythonPackageshome/cplex
~~~

Both of those commands (default and home-specified) invoke the Python package **distutils**. For other options available with that package, consult the documentation of Python **distutils**. 
