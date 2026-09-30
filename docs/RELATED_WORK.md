# Related repositories and selection notes

Four public profiles were reviewed for ideas relevant to a materials-engineering summer research portfolio. Repository descriptions and selected READMEs were examined; this is not a full code audit or a verification of the authors' reported results. The linked repositories remain their authors' work.

| Author / repository | Relevance to this portfolio | What was incorporated |
| --- | --- | --- |
| [AbhinavPJ — Crank–Nicolson TDSE](https://github.com/AbhinavPJ/Crank-Nicholson-Method-for-solving-TDSE) | Closest numerical-method connection: finite differences and a tridiagonal system for a physics equation. Its physical problem is quantum wavepacket evolution, not diffusion. | Topic inspiration for an independently implemented CN comparison on the existing diffusion model. No code, notebook, poster, or data copied. |
| [Rajat Golechha — mesh processing](https://github.com/rajatgolechha3319/graphics_a2_col781) | Geometry processing can support later computational-materials work, but this repository is a graphics assignment. | External reading only. GitHub reports GPL-3.0 and a LICENSE file is present; its code has not been imported. |
| [Divyam Awasthy — Paper-Insight](https://github.com/DivyamAwasthy/Paper-Insight) | Scientific literature retrieval and structured extraction; useful as an example of a research-support tool. | External reading only. Its reported model results are not results of this portfolio. |
| [Priyesha710 — ESP32 code](https://github.com/Priyesha710/ESP32-Code) and [plant-watering system](https://github.com/Priyesha710/plant_watering_system) | Potential instrumentation references. The selected roots contain a firmware sketch and a web application, respectively; the available README descriptions are minimal. | External reading only; no claim about demonstrated materials experiments or hardware performance. |

No license file was visible in the inspected roots of the selected AbhinavPJ, DivyamAwasthy, or Priyesha710 projects, and GitHub's repository metadata did not identify a license. They are linked as references, not redistributed. The GPL mesh-processing project was not copied either because it does not directly implement the transport problem studied here.

Most other visible projects concern games, computer architecture, networking, finance, or general ML. They may demonstrate useful programming techniques, but placing them in this portfolio would not establish materials-research experience.

## The extension made here

The new code solves the *existing* diffusion boundary-value problem with Crank–Nicolson and a Thomas tridiagonal solver. Its own tests compare the tridiagonal calculation against NumPy's dense solver and the concentration profiles against the analytical diffusion series. A parameter sweep compares FTCS and CN error at different time steps. It also demonstrates that linear stability does not guarantee physical concentration bounds for a coarse step.

The equations were taken from the standard numerical method described in [MIT's heat and diffusion notes](https://ocw.mit.edu/courses/18-086-mathematical-methods-for-engineers-ii-spring-2006/resources/am54/), not adapted from another student's source. Code and documentation in this project were prepared with AI assistance, as recorded in [development](DEVELOPMENT.md).
