"""Command line for analysis, verification, figures and isolated run preparation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from .campaign import export,calculate,load_campaign,stage
from .verification import compare_results, verify


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description='HfNbTaTiZr compression data and analysis')
    parser.add_argument('--root',type=Path,default=Path.cwd(),help='Project root (default: current directory).')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('analyze');p.add_argument('--output',type=Path,default=Path('results'));p.add_argument('--overwrite',action='store_true')
    p=sub.add_parser('verify');p.add_argument('--results',type=Path,default=Path('results/analysis.json'))
    p=sub.add_parser('figures');p.add_argument('--output',type=Path,default=Path('docs/figures'))
    p=sub.add_parser('stage');p.add_argument('case_id');p.add_argument('--destination',type=Path,required=True);p.add_argument('--potential-dir',type=Path,required=True)
    args=parser.parse_args(argv);root=args.root.resolve()
    def at_root(path):return path if path.is_absolute() else root/path
    try:
        if args.command=='analyze':
            r=export(root,at_root(args.output),overwrite=args.overwrite)
            print(f"Analyzed {len(r['cases'])} curves; {len(r['sensitivity'])} method-sensitivity configurations.")
            print(f"Saved numerical results to {at_root(args.output)}")
        elif args.command=='verify':
            checked = verify(root, at_root(args.results))
            print(f"PASS: data hashes, {len(checked['cases'])} curves, JSON and four CSV exports; {len(checked['sensitivity'])} processing evaluations.")
        elif args.command=='figures':
            from .plotting import make_figures
            make_figures(root,at_root(args.output));print(f'Figures written to {at_root(args.output)}')
        elif args.command=='stage':
            r=stage(root,args.case_id,at_root(args.destination),at_root(args.potential_dir))
            print(json.dumps(r,ensure_ascii=False,indent=2))
        return 0
    except (ValueError,OSError,KeyError,TypeError,ImportError) as e:
        print(f'ERROR: {e}',file=sys.stderr);return 2
