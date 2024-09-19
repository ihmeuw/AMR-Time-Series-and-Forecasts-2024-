import argparse
import getpass
import pandas as pd
from mcod_prep.utils.mcod_cluster_tools import submit_mcod


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Launcher for map_spec_drug_bug.py")
    parser.add_argument("source", nargs="+")
    write_parser = parser.add_mutually_exclusive_group(required=True)
    write_parser.add_argument(
        '--write', dest='write', action='store_true',
        help='Write the output file from mapping '
        '(only use once your cleaning script has been approved)'
    )
    write_parser.add_argument(
        '--no-write', dest='write', action='store_false',
        help='Do not write the output file from mapping (used for testing purposes)'
    )
    args = parser.parse_args()
    assert {'new'}.isdisjoint(set(args.source)),\
        "The 'new' option doesn't make sense in the launcher yet"
    if 'all' in args.source:
        print("You specified all sources, reading them from nid_metadata...")
        sources = pd.read_csv("FILEPATH").source.unique().tolist()
        active_cols = [c for c in sources if 'active' in c]
        sources = sources.loc[~(sources[active_cols] == 0).all(axis=1)]
    else:
        sources = args.source
    print(f"Submitting jobs for sources {sources}")

    mem = {'SOURCE': "375G", "SOURCE": "50G", "SOURCE": "35G", 'SOURCE': '15G', 'SOURCE': '35G'}
    runtime = {'SOURCE': "05:00:00", "SOURCE": "01:30:00"}
    cores = {'SOURCE': 50}

    worker = "FILEPATH"
    for source in sources:
        jobname = f"map_amr_{source}"
        submit_mcod(
            jobname, language='python', worker=worker,
            cores=cores.get(source, 1), memory=mem.get(source, "10G"),
            runtime=runtime.get(source, '00:30:00'),
            params=[source, '--write'], logging=True
        )
