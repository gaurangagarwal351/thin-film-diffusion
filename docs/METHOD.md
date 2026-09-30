# Model and numerical method

## Physical question

How does the uptake of a diffusing species change with film thickness and diffusivity? This model represents a planar film with lateral dimensions much larger than its thickness. Both faces are exposed to reservoirs that hold the surface concentration at the same constant value.

The scope is idealized transport in a homogeneous material. It is relevant as a starting point for reasoning about diffusion time scales; it does not represent a measured alloy, coating, membrane, or device. A supported film with an impermeable substrate would require a different boundary condition.

![Slab model and numerical stencil](../results/model_schematic.png)

The schematic uses fewer grid nodes for clarity; the saved concentration study uses 201 nodes.

## Equation and assumptions

For constant diffusivity D, Fick's second law is

$$\frac{\partial C}{\partial t}=D\frac{\partial^2 C}{\partial x^2},\qquad 0<x<L.$$

The initial interior concentration is zero. For positive time, C(0,t)=C(L,t)=Cs. Here **L is the full thickness**, not the half-thickness. The reservoirs supply material, so the total species content of the film increases; conservation does not mean constant mass inside this open domain.

With z=x/L, u=C/Cs, and Fo=Dt/L², the equation becomes u_Fo=u_zz on [0,1], with both end values fixed to one. The dimensionless solution is independent of the chosen numerical values of L and D. The diffusion-equation background is covered in [MIT's lecture on Fick's law, conditions, and scaling](https://ocw.mit.edu/courses/3-21-kinetic-processes-in-materials-spring-2006/resources/ls6/).

Assumptions are a constant temperature and diffusivity, no reaction or advection, no swelling, no grain-boundary pathways, and negligible resistance to transfer at either surface. The film does not deplete either external reservoir.

## Finite differences

The solver uses a uniform grid and forward Euler time stepping with a centered second derivative:

$$u_i^{n+1}=u_i^n+r\left(u_{i-1}^n-2u_i^n+u_{i+1}^n\right),\qquad r=\frac{\Delta Fo}{(\Delta z)^2}.$$

For 0<r≤1/2, the update is a weighted combination with non-negative weights. This provides the explicit scheme's stability restriction and preserves the bounds of the initial and boundary values. The study uses r≤0.45 and shortens steps to reach requested output times. See [Strang's MIT notes on the heat and diffusion equation](https://ocw.mit.edu/courses/18-086-mathematical-methods-for-engineers-ii-spring-2006/resources/am54/) for the numerical-method background.

The interior starts at zero and the endpoints at one, representing the imposed surface jump at time 0+. Trapezoidal integration of those nodal values gives an artificial initial mean of 1/(N−1), although the exact continuum initial mean is zero. Saved uptake comparisons therefore begin at Fo=0.002; refinement reduces this initialization effect. Near time zero, the steep boundary layers require particular care.

## Independent analytical calculation

Set v=1−u. It has zero boundary values and initial interior value one. Expanding that constant in a sine series gives coefficients 4/(nπ) for odd n and zero for even n. Each mode decays as exp(−n²π²Fo), so

$$u(z,Fo)=1-\frac{4}{\pi}\sum_{m=0}^{\infty}\frac{\sin[(2m+1)\pi z]}{2m+1}e^{-(2m+1)^2\pi^2 Fo}.$$

Integrating over the thickness gives the mean uptake:

$$\overline{u}(Fo)=1-\frac{8}{\pi^2}\sum_{m=0}^{\infty}\frac{e^{-(2m+1)^2\pi^2 Fo}}{(2m+1)^2}.$$

The implementation evaluates 512 odd modes. At the saved positive times the exponential factors make the omitted tail negligible compared with the discretization error. At Fo=0, the initial condition is returned directly. Extremely early positive times may need more terms; the uptake-time helper only accepts fractions from 1% to below 100%.

## Verification and interpretation

The study compares profiles with the series and refines the grid at Fo=0.05. Because the time step scales as Δz², both the spatial error and the first-order temporal error decrease approximately as Δz². The reported convergence order is therefore for this coupled refinement, not a claim of second-order time integration.

For a chosen uptake, the dimensionless crossing time is fixed. Converting it back using t=Fo L²/D shows the quadratic thickness dependence and inverse diffusivity dependence. The [results report](../results/REPORT.md) states the computed values and errors.

Model agreement is not experimental validation. Applying this calculation to a real film would require material-specific diffusivity, temperature, geometry, surface conditions, and independent measurements. It would also require checking whether the simplifying assumptions are suitable.
