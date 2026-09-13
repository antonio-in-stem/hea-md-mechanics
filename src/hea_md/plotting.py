"""Build figures from the same checked configuration used for numerical extraction."""
from __future__ import annotations
from pathlib import Path
import numpy as np

from .campaign import load_campaign
from .io import load_curve, Curve
from .mechanics import centered_mean
from .verification import verify


def offset_series(curve: Curve, result: dict) -> dict:
    """Expose the plotted offset construction for both figures and tests."""
    cfg = result['config']
    mask = curve.strain <= cfg['analysis_max'] + 1e-12
    window = Curve(curve.strain[mask], curve.stress[mask])
    x, y = centered_mean(window, cfg['smooth_points'])
    e = result['elastic_fit']['slope_GPa']
    b = result['elastic_fit']['intercept_GPa']
    return {'x': x, 'mean': y, 'elastic': e*x+b,
            'offset': e*(x-cfg['offset'])+b, 'raw_x': window.strain,
            'raw_y': window.stress, 'point': result['proof_stress']['point']}


def make_figures(root: Path, output: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root, output = Path(root).resolve(), Path(output).resolve()
    for protected in ('data','src','tests','configs','simulation'):
        if output.is_relative_to(root/protected):
            raise ValueError('Figure output may not overwrite inputs or source.')
    obj = verify(root)
    output.mkdir(parents=True,exist_ok=True)
    cases = load_campaign(root)['cases']
    by_id = {c['id']:c for c in cases}
    curves = {c['id']:load_curve(root/'data/raw'/c['curve_file']) for c in cases}
    results = {r['case_id']:r for r in obj['cases']}
    cfg = obj['analysis_settings']
    maximum = cfg['analysis_max']
    fit_label = f"{100*cfg['fit_min']:g}–{100*cfg['fit_max']:g}%"
    offset_label = f"{100*cfg['offset']:g}%"
    common_label = f"0–{100*maximum:g}%"

    def figure():
        fig = plt.figure(figsize=(7.8,5.1),layout='constrained')
        axes = fig.add_subplot(111)
        axes.tick_params(labelsize=10)
        return fig, axes

    def finish(fig, axes, name, title, xlabel='Engineering compressive strain (%)', ylabel='Engineering compressive stress (GPa)'):
        axes.set_title(title,fontsize=13,pad=12)
        axes.set_xlabel(xlabel,fontsize=11)
        axes.set_ylabel(ylabel,fontsize=11)
        axes.grid(True,alpha=.2)
        axes.spines[['top','right']].set_visible(False)
        with matplotlib.rc_context({'svg.hashsalt':'hfnbtatizr-md'}):
            fig.savefig(output/f'{name}.pdf',metadata={'CreationDate':None,'ModDate':None,'Title':title,
                        'Author':'J. A. Martínez Torres; S. J. Mejía-Rosales'})
            fig.savefig(output/f'{name}.svg',metadata={'Date':None,'Title':title})
            fig.savefig(output/f'{name}.png',dpi=180)
        plt.close(fig)

    fig, ax = figure()
    for c in cases:
        cv = curves[c['id']]
        m = cv.strain <= maximum + 1e-12
        ax.plot(100*cv.strain[m],cv.stress[m],lw=1,label=c['label'])
    ax.set_xlim(0,100*maximum)
    ax.legend(ncols=3,fontsize=9,loc='upper right',frameon=False)
    finish(fig,ax,'compression-common',f'Compression: {common_label} comparison')
    fig, ax = figure()
    for c in cases:
        cv = curves[c['id']]
        ax.plot(100*cv.strain,cv.stress,lw=1,label=c['label'])
    ax.set_xlim(0,100*max(c.strain[-1] for c in curves.values()))
    ax.legend(ncols=3,fontsize=9,loc='upper left',frameon=False)
    finish(fig,ax,'compression-full','Complete compression records')
    for family, ids in [('zr',['eq','zr24','zr28']),('hf',['eq','hf16','hf12']),
                        ('ta',['eq','ta16','ta12']),('nb',['eq','nb24','nb28'])]:
        fig, ax = figure()
        for cid in ids:
            cv = curves[cid]
            m = cv.strain <= maximum + 1e-12
            ax.plot(100*cv.strain[m],cv.stress[m],lw=1.1,label=by_id[cid]['label'])
        ax.set_xlim(0,100*maximum)
        ax.legend(frameon=False,fontsize=10)
        finish(fig,ax,f'composition-{family}',f'{family.capitalize()} composition path')
    fig, ax = figure()
    r = results['eq']
    series = offset_series(curves['eq'],r)
    ax.plot(100*series['raw_x'],series['raw_y'],lw=.8,alpha=.45,label='Recorded stress')
    ax.plot(100*series['x'],series['mean'],lw=1.8,label=f"Centered mean ({cfg['smooth_points']} samples)")
    proof = series['point']
    display_max = min(maximum, max(2.5*cfg['fit_max'],1.5*(proof['strain'] if proof else cfg['fit_max'])))
    line_mask = series['x'] <= display_max
    ax.plot(100*series['x'][line_mask],series['elastic'][line_mask],'--',lw=1.2,label=f'OLS fit, {fit_label}')
    ax.plot(100*series['x'][line_mask],series['offset'][line_mask],':',lw=1.5,label=f'{offset_label} offset line')
    if proof is not None:
        ax.scatter([100*proof['strain']],[proof['stress_GPa']],s=45,zorder=5,label='Operational offset stress')
    ax.set_xlim(0,100*display_max)
    visible_y = series['raw_y'][series['raw_x']<=display_max]
    span = max(float(np.ptp(visible_y)),1e-6)
    ax.set_ylim(float(visible_y.min())-.08*span,float(visible_y.max())+.16*span)
    ax.legend(fontsize=9,frameon=False,loc='upper left')
    finish(fig,ax,'offset-reference','Equiatomic offset construction')
    fig, ax = figure()
    for c in cases:
        rows = [row for row in obj['sensitivity'] if row['case_id']==c['id']
                and row['smooth_points']==cfg['sensitivity_smooth_points'][0]
                and row['confirm_points']==cfg['sensitivity_confirm_points'][0]]
        rows.sort(key=lambda row:row['fit_max'])
        ax.plot([row['fit_max']*100 for row in rows],[row['slope_GPa'] for row in rows],marker='o',ms=3,lw=1.2,label=c['label'])
    ax.margins(y=.28)
    ax.legend(ncols=3,fontsize=9,frameon=False,loc='upper right')
    finish(fig,ax,'stiffness-sensitivity','Fit-window sensitivity',f"Upper fit limit; lower limit {100*cfg['fit_min']:g}% (%)",'Finite-window slope (GPa)')
    fig, ax = figure()
    rows = obj['summary']
    for key, text, marker in [('slope_change_percent',f'{fit_label} stiffness','o'),
                              ('proof_change_percent',f'Operational {offset_label} offset stress','s'),
                              ('peak_change_percent',f'Peak stress, {common_label}','^')]:
        selected = [(i,row[key]) for i,row in enumerate(rows) if row[key] is not None]
        ax.plot([v for _,v in selected],[i for i,_ in selected],marker=marker,ms=5,lw=0,label=text)
    ax.set_yticks(np.arange(len(rows)),[row['label'] for row in rows])
    ax.invert_yaxis()
    ax.axvline(0,ls=':',lw=.7)
    ax.legend(fontsize=9,frameon=False,loc='lower right')
    finish(fig,ax,'relative-properties','Differences from the equiatomic record','Relative change (%)','Nominal composition')
