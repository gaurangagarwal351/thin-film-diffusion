"""Regenerate the documented simulation outputs. Run from the repo root."""

import argparse
import csv
import json
import os
from pathlib import Path
import platform

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from diffusion import analytical_profile, analytical_uptake, physical_time, solve_slab, time_for_uptake


def write_csv(path, headings, rows):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(headings)
        writer.writerows(rows)


def save_figure(fig, directory, name):
    fig.savefig(directory / (name + '.png'), dpi=180, facecolor='white')
    svg_path = directory / (name + '.svg')
    fig.savefig(svg_path, metadata={'Date': None}, facecolor='white')
    svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines()) + '\n')
    plt.close(fig)


def schematic(directory):
    fig, ax = plt.subplots(figsize=(11, 5.3))
    ax.set(xlim=(0, 11), ylim=(0, 5.3))
    ax.axis('off')
    ax.text(0.3, 4.95, 'Model geometry and numerical grid', fontsize=18, weight='bold')
    ax.text(0.3, 4.55, 'Idealized slab exposed to equal concentration reservoirs on both faces', fontsize=11)
    ax.add_patch(Rectangle((3, 2.8), 5, 1.25, facecolor='#e5eff6', edgecolor='#245777', lw=1.7))
    for start, end in [(1.5, 3), (9.5, 8)]:
        ax.annotate('', xy=(end, 3.42), xytext=(start, 3.42), arrowprops={'arrowstyle': '-|>', 'color': '#087e8b', 'lw': 2.2})
    ax.text(0.45, 3.85, 'Reservoir', fontsize=12)
    ax.text(0.45, 3.02, r'$C(0,t)=C_s$', fontsize=13)
    ax.text(8.55, 3.85, 'Reservoir', fontsize=12)
    ax.text(8.55, 3.02, r'$C(L,t)=C_s$', fontsize=13)
    ax.text(5.5, 3.55, 'Uniform, constant diffusivity D', ha='center', fontsize=12)
    ax.text(5.5, 3.1, 'Initial interior concentration: C = 0', ha='center', fontsize=11)
    ax.annotate('', xy=(8, 2.45), xytext=(3, 2.45), arrowprops={'arrowstyle': '<->', 'color': '#245777'})
    ax.text(5.5, 2.1, 'Full film thickness L', ha='center', fontsize=11)
    ax.text(3, 2.13, 'x = 0', ha='center', fontsize=10)
    ax.text(8, 2.13, 'x = L', ha='center', fontsize=10)
    mesh = np.linspace(3, 8, 11)
    ax.plot(mesh, np.ones(11)*1.6, '-o', color='#245777', ms=5)
    ax.scatter([3, 8], [1.6, 1.6], color='#cc6b36', s=55, zorder=3)
    ax.text(0.45, 1.6, 'Uniform grid', va='center', fontsize=11)
    for index, label in [(4, 'i − 1'), (5, 'i'), (6, 'i + 1')]:
        ax.text(mesh[index], 1.18, label, ha='center', fontsize=10)
    ax.text(5.5, 0.64, r'$u_i^{n+1}=u_i^n+r(u_{i-1}^n-2u_i^n+u_{i+1}^n)$', ha='center', fontsize=15)
    ax.text(5.5, 0.16, 'Orange nodes: fixed boundaries   |   Blue nodes: FTCS updates   |   r ≤ 1/2', ha='center', fontsize=10)
    fig.tight_layout(pad=0.5)
    save_figure(fig, directory, 'model_schematic')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results'))
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.alpha': 0.2, 'svg.hashsalt': 'diffusion-study'})

    profile_times = [0.002, 0.01, 0.05, 0.2]
    profiles = solve_slab(profile_times, nodes=201)
    uptake_times = np.linspace(0.002, 0.5, 200)
    uptake_states = solve_slab(uptake_times, nodes=201)
    convergence = []
    for nodes in [21, 41, 81, 161]:
        state = solve_slab([0.05], nodes=nodes)[0]
        error = float(np.max(np.abs(state.concentration - analytical_profile(state.position, state.tau))))
        convergence.append((nodes, nodes-1, 1/(nodes-1), error))
    orders = [float(np.log2(a[3]/b[3])) for a, b in zip(convergence, convergence[1:])]
    tau90 = time_for_uptake(0.9)
    diffusivity = 1e-14  # Illustrative scale only; not a measured material property.
    thicknesses_um = [1, 2, 5]
    t90_seconds = {str(value): float(physical_time(tau90, value*1e-6, diffusivity)) for value in thicknesses_um}
    profile_errors = [float(np.max(np.abs(s.concentration - analytical_profile(s.position, s.tau)))) for s in profiles]
    summary = {
        'data_type': 'simulation; no experimental observations',
        'model': '1D constant-D slab; both surfaces held at normalized concentration 1',
        'nodes': 201, 'maximum_ftcs_ratio': 0.45, 'profile_fourier_times': profile_times,
        'series_terms': 512, 'max_profile_errors': profile_errors,
        'convergence_fourier_time': 0.05, 'observed_orders': orders,
        'fourier_time_to_90_percent': tau90,
        'illustrative_diffusivity_m2_s': diffusivity,
        'time_to_90_percent_seconds_by_thickness_um': t90_seconds,
        'python_version': platform.python_version(), 'numpy_version': np.__version__,
        'matplotlib_version': matplotlib.__version__,
    }
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    write_csv(out/'profiles.csv', ['fourier_time', 'x_over_L', 'C_over_Cs_numerical', 'C_over_Cs_series'],
              ((s.tau, x, y, a) for s in profiles for x, y, a in zip(s.position, s.concentration, analytical_profile(s.position, s.tau))))
    write_csv(out/'uptake.csv', ['fourier_time', 'mean_C_over_Cs_numerical', 'mean_C_over_Cs_series'],
              ((s.tau, s.uptake, analytical_uptake(s.tau)) for s in uptake_states))
    write_csv(out/'grid_convergence.csv', ['nodes', 'intervals', 'normalized_spacing', 'max_abs_error'], convergence)

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.3), layout='constrained')
    colors = ['#156b8a', '#cd7037', '#43865a', '#7c5e9d']
    ax = axes[0, 0]
    for state, color in zip(profiles, colors):
        ax.plot(state.position, state.concentration, color=color, lw=2, label=f'Fo = {state.tau:g}')
        ax.plot(state.position, analytical_profile(state.position, state.tau), '--', color='#222222', lw=0.8, alpha=0.65)
    ax.set(title='Concentration through the film', xlabel='Position x/L', ylabel='Concentration C/Cs', xlim=(0, 1), ylim=(0, 1.04))
    ax.legend(fontsize=9, loc='lower center', ncol=2)
    ax.text(0.5, 0.95, 'Dashed: analytical series', transform=ax.transAxes, ha='center', fontsize=9)

    ax = axes[0, 1]
    ax.plot(uptake_times, [s.uptake for s in uptake_states], color=colors[0], label='Finite differences')
    ax.plot(uptake_times, [analytical_uptake(t) for t in uptake_times], '--', color=colors[1], label='Analytical series')
    ax.axhline(0.9, lw=1, color='#777777', linestyle=':')
    ax.scatter([tau90], [0.9], color=colors[1], zorder=3)
    ax.annotate(f'90% uptake: Fo = {tau90:.4f}', xy=(tau90, 0.9), xytext=(0.15, 0.65), arrowprops={'arrowstyle': '->'})
    ax.set(title='Average uptake', xlabel='Fourier time Fo = Dt/L²', ylabel='Mean C/Cs', xlim=(0, 0.5), ylim=(0, 1.04))
    ax.legend(loc='lower right', fontsize=9)

    ax = axes[1, 0]
    spacing = np.array([row[2] for row in convergence])
    errors = np.array([row[3] for row in convergence])
    ax.loglog(spacing, errors, 'o-', color=colors[0], label='Maximum profile error')
    ax.loglog(spacing, errors[0]*(spacing/spacing[0])**2, '--', color=colors[1], label='Second-order reference')
    ax.set(title='Grid refinement at Fo = 0.05', xlabel='Normalized spacing Δx/L', ylabel='Maximum absolute error in C/Cs')
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    seconds = np.linspace(0, 650, 250)
    for thickness_um, color in zip(thicknesses_um, colors):
        times = diffusivity*seconds/(thickness_um*1e-6)**2
        ax.plot(seconds, [analytical_uptake(t) for t in times], color=color, label=f'L = {thickness_um} µm')
    ax.set(title='Illustrative thickness comparison', xlabel='Time (s)', ylabel='Mean C/Cs', xlim=(0, 650), ylim=(0, 1.04))
    ax.text(0.5, 0.1, 'Assumed D = 10⁻¹⁴ m²/s; not a measured property', transform=ax.transAxes, ha='center', fontsize=9)
    ax.legend(fontsize=9)
    fig.suptitle('Diffusion into a thin film — simulation results', fontsize=16, weight='bold')
    save_figure(fig, out, 'study')
    schematic(out)

    report = f'''# Computed results

Generated by `python run_study.py`. Every number and plot on this page comes from the stated model, not an experiment.

![Simulation results](study.png)

## Numerical checks

The 201-node FTCS solution is compared with a 512-term analytical series. The largest absolute concentration error over the four saved profiles is **{max(profile_errors):.3e}** in normalized concentration units. This covers Fo = 0.002, 0.01, 0.05, and 0.2, not arbitrary times or boundary conditions.

At Fo = 0.05, doubling the number of intervals produces observed convergence orders of **{', '.join(f'{p:.3f}' for p in orders)}**. The time step decreases with the square of grid spacing, so the first-order time error also decreases quadratically with spacing. This check measures discretization accuracy against the same mathematical model; it does not validate a real material.

## Thickness study

The analytical solution reaches 90% mean uptake at **Fo = {tau90:.6f}**. With an illustrative D = 1e-14 m²/s:

| Film thickness | Time to 90% uptake |
| --- | --- |
''' + '\n'.join(f'| {thickness} µm | {t90_seconds[str(thickness)]:.3f} s |' for thickness in thicknesses_um) + '''

At fixed diffusivity, doubling thickness multiplies the time to a given uptake by four because t = Fo L²/D. This is a consequence of the model's scaling, not an experimentally discovered relation. No specific film composition or temperature is assigned to the illustrative diffusivity.

## Files and interpretation

- `profiles.csv`: position, numerical concentration, and series concentration at four times.
- `uptake.csv`: numerical and analytical mean uptake over time.
- `grid_convergence.csv`: maximum error for 21, 41, 81, and 161 nodes.
- `summary.json`: parameters, numerical checks, and software versions.
- `model_schematic.svg`: idealized slab geometry and finite-difference stencil.

There are no repeats, uncertainty bars, or fitted material parameters because these outputs are deterministic model calculations. The model omits reactions, microstructure, concentration-dependent diffusivity, and finite surface-transfer resistance. See the [method](../docs/METHOD.md) and [development note](../docs/DEVELOPMENT.md).
'''
    (out/'REPORT.md').write_text(report)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
