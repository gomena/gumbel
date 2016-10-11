---
layout: post
title: "Code for On quadrature methods for refractory point processes"
categories: code
author: Gonzalo Mena
excerpt_separator: <!--more-->
comments: true
---


[In the GitHub site](https://github.com/gomena/RefractoryLikelihoods) there are two examples (in Matlab) for illustrating the methods presented in the article.
The first example (ExampleRenewal.m) corresponds to a Renewal Process (with inverse Gaussian ISI-distribution and an absolute refractory period). The log likelihood is computed using the four presented methods (DR1, DR2, CT, GL)
The second example (ExampleGLM.m) is a simple GLM model with one dimensional covariate, two parameters to be estimated and a piecewise linear refractory function r(t). First, a spike train is simulated (T=40s), and then estimation is done using Newton Method with backtracking linesearches, using both DR1 and GL.
Nodes and weights were already saved in .mat files. There are several ways of computing these weights. For the GL quadrature, We used [Greg von Winckel's algorithm](http://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights/content/lglnodes.m) , in which the nodes are obtained as the roots of the derivatives of the Legendre polynomials, and the weights are a function of the evaluations of the polynomials at the corresponding nodes
See *C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods and in Fluid Dynamics," Section 2.3. Springer-Verlag 1987* for details


[These are the slides](These are the slides) of an oral presentation of this work at the [2014 Minghui Yu Memorial Conference](http://www.columbia.edu/~bmr2136/MYMC2014/) , organized by the Statistics Department, Columbia University.
