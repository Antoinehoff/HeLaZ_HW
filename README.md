This is a solver of the Hasegawa-Wakatani equations :

<img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;\partial_t&space;\zeta&space;&plus;&space;\{\phi,\zeta\}&space;&=&space;\alpha&space;(\phi&space;-&space;n)&space;-&space;\mu\nabla^4\zeta\nonumber\\&space;\partial_t&space;n&space;&plus;&space;\{\phi,n\}&space;&=&space;\alpha&space;(\phi&space;-&space;n)&space;-&space;\mu\nabla^4&space;n&space;-\kappa\partial_y\phi&space;\end{align*}" title="\begin{align*} \partial_t \zeta + \{\phi,\zeta\} &= \alpha (\phi - n) - \mu\nabla^4\zeta\nonumber\\ \partial_t n + \{\phi,n\} &= \alpha (\phi - n) - \mu\nabla^4 n -\kappa\partial_y\phi \end{align*}" />

using a nonlinear Fortran spectral solver. The structure of the code is highly inspired by GBS and the results have been compared with the work of Ammar Hakim http://ammar-hakim.org/sj/je/je17/je17-hasegawa-wakatani.html

To run it :
Launch matlab, run setup.m, then run.m then analysis.m. You can also launch the code manually with ./../bin/helaz from wk/ with an appropriate set of parameters put in wk/fort.90


![](zeta_00.gif)
