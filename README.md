# Heat equation in 2D
Solved using semi-implicit time discretization. The code allows inclusion of any type and number of "holes" in the rectangular region
and it is assumed that no diffusion occurs across holes (they're insulated). 

The left hand side of the equation gives the change of the temperature of the point in time.

The right hand side is the diffusivity times the Lapacian of the temperature in space. The Laplacian is the divergence of the gradiaent. So we
get the divergence of the gradient (in space) of temperature. Divergence of a vector field is the extent to which the vector field flux behaves like a source at a given point. 

The code files are organised as follows:
1. **heat.m** contains code to implement heat diffusion on a rectangular grid. 
2. **dirichlet.m** conatins code for the Dirichlet bcs. **neumann.m** similarly. 
3. **location.m** contains code to look for closest point inside the region. 

The gif files are configured to 2 holes (1 circle, 1 rectangle) and 3 initial sources of different temperatures for Dirichlet and Neumann boundaries.  

Very useful pages:
<ul>
  <li>https://en.wikipedia.org/wiki/Divergence</li>
  <li>https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf</li>
</ul>



Thanks for visiting :)
 
