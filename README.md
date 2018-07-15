# Repeat aware scaffolding evaluation framework by Igor Mandric


### General information

Genome scaffolding is a classical challenging problem in bioinformatics. It refers to joining assembly contigs into chains (called scaffolds). The join between two contigs A and B is considered correct if:

  * Their relative orientation is correct
  * Their relative order is correct
  * The gap estimate is similar to the true distance on the reference
  
The problem of scaffolding validation is also a challenging one. One of the main issues which hinders from an adequate scaffolding evaluation are genome repeats. The previous standard for evaluation [(Hunt et al., *Genome Biology*, 2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r42) did not take into account repeats. In this evaluation framework, repeats are taken into account. 
  
<p align="center">
  <img src="http://alan.cs.gsu.edu/repeat-aware/figure.png" width="75%">
</p>

The new evaluation framework considers the optimal assignment of contigs in the output scaffolding to contigs in the reference scaffolding in the sense of the number of correct links.



### Software prerequisites

Our framework is written in [**Python**](https://www.python.org/) programming language. By default, **Python** is available on any Linux platform. The optimal assignment of **Inferred Scaffolding** to **Reference Scafolding** is performed by using [**IBM CPLEX**](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) optimizer.

We use the following scripts/libraries:

- [Biopython](http://biopython.org/) library
- Python API of CPLEX
- [MUMmer](http://mummer.sourceforge.net/) (rapid alignment of very large DNA and amino acid sequences)

##### Install Biopython

To install Biopython, type the following commands in your prompt:

~~~bash
pip install numpy
pip install biopython
~~~

##### Install Python API of CPLEX

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

##### Install MUMmer

The current **MUMmer** release can be downloaded from [SourceForge](http://sourceforge.net/projects/mummer). To install **MUMmer**, type the following commands:

~~~bash
tar -xvzf MUMmer3.0.tar.gz
make check
make install
~~~

Our framework uses **nucmer**, **delta-filter**, and **show-coords** scripts from **MUMmer** package.



### Running the framework

Our framework takes as input the following:

- Reference genome fasta file
- Assembly contigs i.e., contigs obtained by running an assembly program. For example:
  - Velvet
  - ALLPATHS-LG
  - SOAPdenovo2
- Scaffolding output fasta file i.e. scaffolds obtained by running a scaffolding tool. For example:
  - ScaffMatch
  - OPERA-LG
  - BESST

We first have to prepare **Reference Scaffolding** dataset which is based on real assembly contigs. Correction of misassemblies has to be done and the output scaffolding is written to the **ref_scaf** file in the **$output** directory:

~~~bash
python build_ref_scaf.py --reference $reference --contigs $assembly_contigs --prefix "$outdir/$ref_scaf"
~~~

**ref_scaf.scaf** is a file of **.scaf** format similar to one used by [**OPERA**](https://sourceforge.net/projects/operasf/) tool.

Run your tool of choice with the contigs from the file "$outdir/$ref_scaf" and get the scaffolds. 
Now, having the output scaffolding, we have to create the **out_scaf** **.scaf** file which corresponds to the **Inferred Scaffolding**. In our framework, we align the assembly contigs to the output scaffolding using **nucmer**. Run the following command to get the **out_scaf** file:

~~~bash
python build_out_scaf.py --scafolds $scaffolds --contigs "$outdir/ref_scaf.fa" --outscaf "$outdir/out_scaf"
~~~

Finally, to evaluate a scaffolding tool we have to find the assignment of the **Inferred Scaffolding** contigs to the **Reference Scaffolding** contigs. The following command peforms the optimal assignment i.e the assignment which keeps the maximum number of contig links:

~~~bash
python validation.py --outscaf "$outdir/out_scaf" --refscaf "$outdir/ref_scaf"
~~~

The output is the number of correct contig links and other measures - sensitivity and PPV (positive predictive value).


### Help

  * In order to get more help for the validation.py script, execute the command:
   ```
  python validation.py -h
  ```
   You should see the following output:
   ```
usage: validation.py [-h] --outscaf OUTSCAF --refscaf REFSCAF
                            [--dist DISTANCE_THRESHOLD] [--noskip NOSKIP]

Repeat-aware evaluation framework of Igor Mandric

optional arguments:
  -h, --help            show this help message and exit
  --outscaf OUTSCAF     output scaffolding
  --refscaf REFSCAF     reference scaffolding
  --dist DISTANCE_THRESHOLD
                        links for which estimated gap differs from the true
                        value by more than DISTANCE_THRESHOLD bp are
                        considered wrong (default 1000 bp)
  --noskip NOSKIP       if the argument is present than do not consider links
                        skipping contigs as wrong if the contigs are at most
                        NOSKIP bp apart on the reference, otherwise links
                        skipping contigs are wrong

   ```



### Contact us

In order to get more information on the framework, do not hesitate to contact [Igor Mandric (Georgia State University)](mailto:mandric.igor@gmail.com) or visit the website [link](https://mandricigor.github.io/repeat-aware/).
