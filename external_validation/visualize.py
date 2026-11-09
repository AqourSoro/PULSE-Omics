"""Publication-quality Seaborn plots for independent PULSE-Omics validation.

Reads previously generated real experimental results; never retrains models,
alters raw data, or fabricates an unreported paper benchmark statistic.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DEFAULT_ROOT = Path('/mnt/f/datasets/pulse-omics')
PAPER_DOI = '10.1038/s43588-026-01026-5'
PAPER_REFERENCE = Path(__file__).with_name('paper_reference.csv')
ARM_NAMES = {'pulse_full': 'PULSE + history',
             'no_history_input': 'PULSE − history',
             'no_alignment': 'PULSE − alignment',
             'pulse_current': 'PULSE current-only',
             'last_observation': 'Last observation'}
PALETTE = {'PULSE + history': '#216e99',
           'PULSE − history': '#f0a64a',
           'PULSE − alignment': '#52a88d',
           'PULSE current-only': '#9b80b8',
           'Last observation': '#a0a8b0'}
PAPER_ARMS = {'pulse_full': 'PULSE (baseline + follow-up labs)',
              'past_only': 'PULSE (baseline full profile)',
              'followup_only': 'PULSE (follow-up labs only)',
              'midas': 'MIDAS (follow-up labs only)',
              'scvaeit': 'scVAEIT (follow-up labs only)',
              'stabmap': 'StabMap (follow-up labs only)'}


def latest_run(root: Path) -> Path:
    base = root/'results'/'pulse_validation'
    found = sorted((p for p in base.iterdir() if p.is_dir() and (p/'summary.csv').is_file()),
                   key=lambda p: p.name) if base.exists() else []
    if not found:
        raise FileNotFoundError(f'No completed summary.csv found in {base}')
    return found[-1]


def style():
    sns.set_theme(style='whitegrid', context='notebook', font_scale=1.05)
    plt.rcParams.update({'figure.dpi': 135, 'savefig.dpi': 220,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'font.family': 'DejaVu Sans',
                         'figure.facecolor': 'white', 'axes.titleweight': 'bold'})


def save(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / (stem + '.png'), bbox_inches='tight', facecolor='white')
    fig.savefig(outdir / (stem + '.svg'), bbox_inches='tight', facecolor='white')
    plt.close(fig)


def load_run(run_dir: Path):
    summary = pd.read_csv(run_dir/'summary.csv')
    required = {'dataset','source','target','task','arm','status','pearson_median','mae_micro'}
    if not required.issubset(summary.columns):
        raise ValueError('Run summary lacks fields: '+str(sorted(required-set(summary.columns))))
    summary = summary.loc[summary.status=='completed'].copy()
    if summary.empty:
        raise ValueError(f'No completed experimental arms in {run_dir}')
    for field in ('pearson_median', 'mae_micro', 'mse_micro', 'change_pearson_median'):
        if field in summary:
            summary[field] = pd.to_numeric(summary[field], errors='coerce')
    summary['direction'] = summary.source.map({'protein':'Protein → Metabolite','metab':'Metabolite → Protein'})
    summary['model'] = summary.arm.map(ARM_NAMES).fillna(summary.arm)
    return summary


def panel_overview(summary, outdir, run_dir):
    directions = list(summary.direction.dropna().unique())
    fig, axes=plt.subplots(len(directions), 2, figsize=(15, max(4.8, 4.2*len(directions))), squeeze=False)
    for row, direction in enumerate(directions):
        frame = summary.loc[summary.direction==direction].copy()
        frame['task_model'] = frame['dataset'].str.upper()+' · '+frame['task'].map({'current_completion':'Current only',
                            'historical_completion':'Historical target'}).fillna(frame['task']) + ' · ' + frame['model']
        frame = frame.sort_values(['task', 'arm'])
        for col, (field, label) in enumerate([('pearson_median', 'Per-feature Pearson r · median ↑'),
                                               ('mae_micro', 'MAE · observed target units ↓')]):
            ax=axes[row, col]
            valids = frame.dropna(subset=[field]).copy()
            if valids.empty:
                ax.text(.5,.5,'No values',transform=ax.transAxes,ha='center')
                continue
            colors=[PALETTE.get(x,'#607c8e') for x in valids.model]
            bars=ax.barh(valids.task_model, valids[field],color=colors,height=.68,alpha=.94)
            pad=max(abs(valids[field]).max()*.025,.008)
            for b,v in zip(bars,valids[field]):
                ax.text(v+pad,b.get_y()+b.get_height()/2, f'{v:.3f}', va='center',fontsize=9)
            ax.set_xlim(min(0., float(valids[field].min())-.05),max(0.,float(valids[field].max())+pad*8))
            ax.set_xlabel(label)
            ax.set_title(f'{direction} | {label.split(" · ")[0]}')
            sns.despine(ax=ax,left=True,bottom=False)
    note = ('Different eligibility/splits in legacy current-only versus historical tasks: '
            'do NOT interpret their difference as history gain. '
            'Full vs ablation within historical task uses paired test samples.')
    fig.suptitle('PULSE independent-cohort validation', fontsize=17, fontweight='bold', y=1.02)
    fig.text(.02,-.015,note,fontsize=9,color='#566472',wrap=True)
    fig.tight_layout(h_pad=3,w_pad=4)
    save(fig,outdir,'01_model_overview')


def read_feature_metrics(run_dir,summary,task,filename):
    chunks=[]
    for _, row in summary.loc[summary.task==task].iterrows():
        path=run_dir/str(row.dataset)/f'{row.source}_to_{row.target}'/str(row.task)/str(row.arm)/filename
        if not path.is_file():
            continue
        frame=pd.read_csv(path)
        if not {'feature','pearson_r','mae','mse'}.issubset(frame.columns):
            raise ValueError(f'Invalid feature metrics columns: {path}')
        frame['dataset']=row.dataset
        frame['direction']=row.direction
        frame['arm']=row.arm
        frame['model']=row.model
        chunks.append(frame)
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def feature_distribution(frame,outdir):
    if frame.empty:
        return
    frame=frame.loc[frame.arm.isin(['pulse_full','no_history_input','no_alignment','last_observation'])].copy()
    frame['group']=frame.dataset.astype(str)+' / '+frame.direction.astype(str)
    groups=frame.group.unique().tolist()
    if not groups:
        return
    fig,axes=plt.subplots(len(groups),1,figsize=(13, max(4.5,4.3*len(groups))),squeeze=False)
    order=['PULSE + history','PULSE − history','PULSE − alignment','Last observation']
    for i,group in enumerate(groups):
        ax=axes[i,0]
        subset=frame.loc[(frame.group==group)&frame.pearson_r.notna()]
        models=[x for x in order if x in subset.model.unique()]
        if not models:
            ax.text(.5,.5,'No valid per-feature correlations',ha='center',va='center',transform=ax.transAxes)
            continue
        sns.boxplot(data=subset,x='model',y='pearson_r',hue='model',order=models,
                    hue_order=models,palette=PALETTE,legend=False,
                    showfliers=False,width=.58,ax=ax)
        sns.stripplot(data=subset,x='model',y='pearson_r',order=models,
                      color='#263547',alpha=.16,size=2.5,jitter=.22,ax=ax)
        ax.axhline(0,color='#94a3b8',lw=.8)
        ax.set_xlabel('')
        ax.set_ylabel('Pearson r per feature')
        ax.set_title(f'{group} · held-out features')
        for label in ax.get_xticklabels():label.set_rotation(7)
        ax.set_ylim(-1.04,1.04)
        sns.despine(ax=ax)
    fig.suptitle('Per-feature accuracy distributions · matched historical task',fontsize=16,fontweight='bold',y=1.01)
    fig.tight_layout(h_pad=3)
    save(fig,outdir,'02_feature_correlations')


def ablation_scatter(frame,outdir):
    if frame.empty:
        return pd.DataFrame()
    subsets=[]
    for (dataset,direction), group in frame.groupby(['dataset','direction'],sort=True):
        for comparison in ('no_history_input','no_alignment'):
            a=group.loc[group.arm=='pulse_full',['feature','pearson_r']].rename(columns={'pearson_r':'pulse_r'})
            b=group.loc[group.arm==comparison,['feature','pearson_r']].rename(columns={'pearson_r':'baseline_r'})
            joined=a.merge(b,on='feature').dropna()
            if joined.empty:
                continue
            joined['dataset']=dataset
            joined['direction']=direction
            joined['comparison']=ARM_NAMES[comparison]
            joined['gain_r']=joined.pulse_r-joined.baseline_r
            subsets.append(joined)
    if not subsets:
        return pd.DataFrame()
    paired=pd.concat(subsets,ignore_index=True)
    paired.to_csv(outdir/'paired_feature_effects.csv',index=False)
    cats=list(paired[['dataset','direction','comparison']].drop_duplicates().itertuples(index=False,name=None))
    ncols=2
    nrows=(len(cats)+1)//2
    fig,axes=plt.subplots(nrows,ncols,figsize=(13,5*nrows),squeeze=False)
    for i,(ds,direction,contrast) in enumerate(cats):
        ax=axes.flat[i]
        view=paired.loc[(paired.dataset==ds)&(paired.direction==direction)&(paired.comparison==contrast)]
        sns.scatterplot(data=view,x='baseline_r',y='pulse_r',alpha=.6,s=30,color='#216e99',edgecolor='none',ax=ax)
        ax.plot([-1,1],[-1,1],linestyle='--',color='#78838e',lw=1)
        ax.set(xlim=(-1.03,1.03),ylim=(-1.03,1.03),aspect='equal',
               xlabel=contrast+' · Pearson r',ylabel='Full PULSE · Pearson r')
        ax.set_title(f'{ds} · {direction}\nMedian feature Δr = {view.gain_r.median():+.3f}   n={len(view)}')
        sns.despine(ax=ax)
    for ax in axes.flat[len(cats):]:ax.axis('off')
    fig.suptitle('Matched-feature ablation comparisons · historical task',fontsize=16,fontweight='bold',y=1.01)
    fig.tight_layout(h_pad=2.5,w_pad=3)
    save(fig,outdir,'03_paired_feature_effects')
    return paired


def paper_comparison(summary,reference,outdir):
    paper=pd.read_csv(reference)
    required={'figure','task','metric','model','value','population','aggregation','source'}
    if not required.issubset(paper.columns):
        raise ValueError('Paper reference file missing columns: '+str(sorted(required-set(paper.columns))))
    fig2=paper.loc[(paper.figure=='Fig. 2c') & (paper.metric=='pearson_r')].copy()
    if fig2.empty:
        raise ValueError('Reference lacks verified Fig. 2c Pearson r values')
    fig2['model_short']=fig2.model
    ours=summary.loc[(summary.task=='historical_completion') & (summary.target=='metab') &
                     (summary.arm.isin(['pulse_full','no_history_input','no_alignment']))].copy()
    ours['group']=ours.dataset.astype(str)+' · '+ours.direction.astype(str)+' · '+ours.model.astype(str)
    if 'change_pearson_median' not in ours:
        ours['change_pearson_median']=np.nan
    ours=ours.dropna(subset=['change_pearson_median'])
    fig,axes=plt.subplots(1,2,figsize=(17,7))
    bars=axes[0].barh(fig2.model_short,fig2.value,color=['#216e99' if x.startswith('PULSE full') else '#95aac1' for x in fig2.model_short],height=.66)
    for b,v in zip(bars,fig2.value): axes[0].text(v+.014,b.get_y()+b.get_height()/2,f'{v:.3f}',va='center')
    axes[0].set_xlim(0,1.05);axes[0].set_xlabel('Pearson r · change in ONE LDL lipid marker')
    axes[0].set_title('Original paper · UK Biobank · Fig. 2c')
    if not ours.empty:
        ours=ours.sort_values(['dataset','direction','arm'])
        sns.barplot(data=ours,y='group',x='change_pearson_median',hue='model',dodge=False,
                    palette=PALETTE,legend=False,ax=axes[1],errorbar=None)
        for bar,v in zip(axes[1].patches,ours.change_pearson_median):
            if np.isfinite(v):axes[1].text(max(v,0)+.015,bar.get_y()+bar.get_height()/2,f'{v:.3f}',va='center',fontsize=9)
    else:axes[1].text(.5,.5,'No valid change correlations',ha='center',va='center',transform=axes[1].transAxes)
    axes[1].set_xlim(-.15,1.1);axes[1].set_xlabel('Median per-feature change Pearson r')
    axes[1].set_title('Independent cohorts · matched historical task')
    fig.suptitle('Same metric family, DIFFERENT biological endpoint and aggregation',fontsize=16,fontweight='bold',y=1.02)
    fig.text(.02,-.015,'Left: UKB change in one LDL lipid biomarker. Right: independent-cohort metabolite change correlations, median over measured features. These are NOT directly equivalent.',
             fontsize=9,color='#566472')
    fig.tight_layout(w_pad=4)
    save(fig,outdir,'04_paper_context')
    rows=[]
    for _,r in fig2.iterrows():
        rows.append({'origin':'PULSE paper','dataset':r.population,'figure':r.figure,'model':r.model,
                     'metric':'Pearson r, longitudinal change','aggregation':r.aggregation,'value':r.value,
                     'source':r.source,'comparability':'Reference only: UKB single LDL marker'})
    for _,r in ours.iterrows():
        rows.append({'origin':'Independent experiment','dataset':r.dataset,'figure':'local experiment',
                     'model':str(r.direction)+' / '+str(r.model),'metric':'Pearson r, longitudinal change',
                     'aggregation':'median across assay features','value':r.change_pearson_median,
                     'source':'run summary.csv','comparability':'Not directly comparable: different biomarkers/assays'})
    table=pd.DataFrame(rows)
    table.to_csv(outdir/'paper_comparison.csv',index=False)
    return table


def write_report(summary,comparison,paired,outdir,run_dir):
    lines=['# Independent external validation · graphical report','',
           f'Experimental run: `{run_dir}`','',
           '## Experimental results','',
           '| Dataset | Direction | Task | Model | Median Pearson r | MAE | MSE | Change median r |',
           '|---|---|---|---|---:|---:|---:|---:|']
    def f(x):
        try:
            v=float(x)
            return f'{v:.4f}' if np.isfinite(v) else 'NA'
        except (ValueError,TypeError):return 'NA'
    for _,r in summary.iterrows():
        lines.append(f'| {r.dataset} | {r.direction} | {r.task} | {r.model} | {f(r.pearson_median)} | '
                     f'{f(r.mae_micro)} | {f(r.get("mse_micro"))} | {f(r.get("change_pearson_median"))} |')
    lines+=['','## Verified paper reference · Fig. 2c','',
            'The paper compares longitudinal changes for a **single free-cholesterol / total-lipid ratio in large LDL particles**, while our values are medians across available assayed features. Direct cross-study accuracy ranking is **not valid**.','',
            '| Reference model | Paper Pearson r |','|---|---:|']
    for _,r in comparison.loc[comparison.origin=='PULSE paper'].iterrows():
        lines.append(f'| {r.model} | {f(r.value)} |')
    lines+=['','Paper DOI: https://doi.org/'+PAPER_DOI,
            'Original paper Fig. 2a reports all 251 metabolites; its precise table-level model medians are not included here because they have not been extracted from the official Supplementary Table 2. They must not be inferred from a plot.',
            'Fig. 2c Pearson r values are transcribed from the annotated printed scatter panels of the original article.',
            '', '## Within-cohort comparisons','']
    if not paired.empty:
        for (ds,direction,contrast), group in paired.groupby(['dataset','direction','comparison']):
            lines.append(f'- {ds}, {direction}, full versus {contrast}: median per-feature Δr = '
                         f'{group.gain_r.median():+.4f} across {len(group)} matched features.')
    else:lines.append('No valid paired feature contrasts.')
    lines+=['','### Interpretation','',
            'Only historical-task full-versus-ablation contrasts are matched on evaluation samples. Old run files may use different eligibility/splits for current-only versus historical-task arms.',
            'Cross-cohort raw MAE/MSE values are assay/unit dependent. Per-feature Pearson r and longitudinal-change correlation are the most interpretable shared metric families, but do not by themselves equate the biological tasks.',
            'No multi-seed uncertainty or significance has been calculated in this single-run report.','']
    (outdir/'comparison_report.md').write_text('\n'.join(lines),encoding='utf-8')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--run-dir',type=Path,help='Existing results/.../<UTC run ID>; default newest')
    p.add_argument('--out-dir',type=Path,help='Default <run-dir>/figures')
    p.add_argument('--paper-reference',type=Path,default=PAPER_REFERENCE)
    args=p.parse_args(argv)
    run_dir=args.run_dir or latest_run(args.data_root)
    outdir=args.out_dir or run_dir/'figures'
    outdir.mkdir(parents=True,exist_ok=True)
    style()
    summary=load_run(run_dir)
    panel_overview(summary,outdir,run_dir)
    feature=read_feature_metrics(run_dir,summary,'historical_completion','feature_metrics.csv')
    feature_distribution(feature,outdir)
    paired=ablation_scatter(feature,outdir)
    comparison=paper_comparison(summary,args.paper_reference,outdir)
    write_report(summary,comparison,paired,outdir,run_dir)
    print('Saved plot and comparison reports to:',outdir)
    for path in sorted(outdir.iterdir()):print(' -',path.name)
    return 0


if __name__=='__main__':
    raise SystemExit(main())