![](zeta_00.gif)

This is a solver of the Hasegawa-Wakatani equations :

<img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;\partial_t&space;\zeta&space;&plus;&space;\{\phi,\zeta\}&space;&=&space;\alpha&space;(\phi&space;-&space;n)&space;-&space;\mu\nabla^4\zeta\nonumber\\&space;\partial_t&space;n&space;&plus;&space;\{\phi,n\}&space;&=&space;\alpha&space;(\phi&space;-&space;n)&space;-&space;\mu\nabla^4&space;n&space;-\kappa\partial_y\phi&space;\end{align*}" title="\begin{align*} \partial_t \zeta + \{\phi,\zeta\} &= \alpha (\phi - n) - \mu\nabla^4\zeta\nonumber\\ \partial_t n + \{\phi,n\} &= \alpha (\phi - n) - \mu\nabla^4 n -\kappa\partial_y\phi \end{align*}" />

using a nonlinear Fortran spectral solver. The structure of the code is highly inspired by GBS and the results have been compared with the work of Ammar Hakim http://ammar-hakim.org/sj/je/je17/je17-hasegawa-wakatani.html

To run it :
  open local/dirs.inc and change the path of the main folder accordingly to yours.

  type _make_ in HeLaz/

  launch matlab to run setup.m (you can easily tweak the parameters from them) that writes automatically the input parameters file wk/fort.90

  run.m to run it or you can also launch the code manually with _./../bin/helaz_ from wk/

  you can then use analysis.m to obtain various plots and gifs

