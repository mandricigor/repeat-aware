
<h1 align="center">
Repeat Aware Evaluation of Scaffolding Tools
</h1>


  [Home](index.md) |
  [Install&Use](install.md) |
  [Datasets](datasets.md) |
  [About](about.md)

<link rel="stylesheet" type="text/css" href="http://alan.cs.gsu.edu/repeat-aware/mystyle.css">


##### Share this page on social networks: 
[<img alt="Share on Facebook" src="http://alan.cs.gsu.edu/repeat-aware/images/flat_web_icon_set/black/Facebook.png">](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fmandricigor.github.io%2Frepeat-aware%2F&t=)
[<img alt="Share on Google+" src="http://alan.cs.gsu.edu/repeat-aware/images/flat_web_icon_set/black/Google+.png">](https://plus.google.com/share?url=https%3A%2F%2Fmandricigor.github.io%2Frepeat-aware%2F)
[<img alt="Submit to Reddit" src="http://alan.cs.gsu.edu/repeat-aware/images/flat_web_icon_set/black/Reddit.png">](http://www.reddit.com/submit?url=https%3A%2F%2Fmandricigor.github.io%2Frepeat-aware%2F&title=)
[<img alt="Share on LinkedIn" src="http://alan.cs.gsu.edu/repeat-aware/images/flat_web_icon_set/black/LinkedIn.png">](http://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fmandricigor.github.io%2Frepeat-aware%2F&title=&summary=&source=https%3A%2F%2Fmandricigor.github.io%2Frepeat-aware%2F)


The aim of this project is to provide a framework for **repeat aware** evaluation of scaffolding tools. The first comprehensive scaffolding evaluation was performed in [(Hunt et al., *Genome Biology*, 2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r42). Its main drawback is that it considers only the **"best match"** for each contig, i.e. the alignment of the contig to the genome with the highest similarity score.

We proposed a new evaluation framework which has the following advantages:

- Repeat-awareness
- Fair evaluation
- Ability to evaluate repeat-aware scaffolders (like [OPERA-LG](https://sourceforge.net/projects/operasf/))
- Ease-of-use


In the picture below (credit to Sergey Knyazev, GSU) is a simple example illustrating the main idea.
  
<p align="center">
  <img src="http://alan.cs.gsu.edu/repeat-aware/figure.png">
</p>

- **Reference Scaffolding** is the "golden true" scaffolding
- **Inferred Scaffolding** is the scaffolding produced by a tool

**Contig 2** has two copies in the reference, namely **2a** and **2b**. **Contig 4** has only one copy, but is inferred to have two copies **4a** and **4b** in the inferred scaffolding.


Our approach is to assign contigs in the **Inferred Scaffolding** to the contigs of the **Reference Scaffolding** maximizing the number of correct links. Asigning **Contig 2** and **Contig 4a** from the **Inferred Scaffolding** to **Contig 2a** and **Contig 4** in the **Reference Scaffolding** correspondingly, we obtaing **2** correct links. Any other assingment delivers less correct links.

