# Understanding and extending the study

This is a short guide for working through the code and explaining its limits.

1. **Run the tests and study script.** Open `results/study.png` and compare it with the committed figure. The numerical curves should closely follow the analytical ones.
2. **Read the boundary conditions first.** Both faces are held at the same concentration. This creates symmetric profiles, with the lowest concentration at the film centre. A film bonded to an impermeable substrate would be a different problem.
3. **Follow one update in `solve_slab`.** Each interior value moves toward the average of its neighbours. The two endpoints remain fixed. The update uses the old values on its entire right-hand side before modifying the interior array.
4. **Explain the stability restriction.** The update weights are r, 1−2r, and r. They stay non-negative for r≤1/2. The function rejects larger values rather than silently producing an unstable calculation.
5. **Distinguish model verification from experimental validation.** The series solution and numerical solver describe the same assumed physical system. Their agreement checks the calculation. It does not prove that the assumptions fit a particular material.
6. **Check the thickness scaling.** The normalized solution depends on Fo=Dt/L². Double L with D fixed: the same uptake takes four times as long. Double D with L fixed: it takes half as long.
7. **Explain the convergence figure.** The grid spacing and time step decrease together. A slope near two is expected for this coupled refinement; the time integrator itself is first order.
8. **Name the largest limitation.** No measured D, film composition, or experiment is used. A real comparison needs those inputs and a reason to assume constant diffusivity and instantaneous surface equilibration.

## Small exercises

- Change the thicknesses in `run_study.py` and predict the 90% uptake times before running it.
- Use `solve_slab` with a finer grid at Fo=0.002 and check how the largest profile error changes.
- Compare `analytical_profile` using 4, 16, and 512 terms at early times. Look for the effect of series truncation near the boundary.

Record the parameters, your prediction, the actual output, and an explanation when you complete an exercise. These exercises are suggestions, not completed work attributed to the repository owner.
