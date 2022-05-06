# Heat equation in 2D
Solved using semi-implicit time discretization. The code allows inclusion of any type and number of "holes" in the rectangular region
and it is assumed that no diffusion occurs across holes (they're insulated). 

$`\sqrt{2}`$.

The left hand side of the equation

The code files are organised as follows:
1. **heat.m** contains code to implement heat diffusion on a rectangular grid. 
2. **dirichlet.m** conatins code for the Dirichlet bcs. **neumann.m** similarly. 
3. **location.m** contains code to look for closest point inside the region. 

Very useful pages:
<ul>
  <li>https://en.wikipedia.org/wiki/Divergence</li>
  <li>https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf</li>
</ul>



Thanks for visiting :)
 
