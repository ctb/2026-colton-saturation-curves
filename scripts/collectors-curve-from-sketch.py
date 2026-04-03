#! /usr/bin/env python
import sys
import csv
import argparse
import sourmash
import random
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sketch')
    p.add_argument('-o', '--output-csv', required=True)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=1000, type=int)
    p.add_argument('-N', '--sampling-interval', type=int, default=10_000)
    args = p.parse_args()

    db = sourmash.load_file_as_index(args.sketch)
    db = db.select(ksize=args.ksize, scaled=args.scaled)
    if len(db) != 1:
        print(f'ERROR: {len(db)} signatures found. Need exactly one.',
              file=sys.stderr)
        sys.exit(-1)

    ss = list(db.signatures())[0]
    mh = ss.minhash.downsample(scaled=args.scaled)
    print(f"Loaded signature '{ss.name}' - len {len(mh)} / sum {mh.sum_abundances}")
    if not mh.track_abundance:
        print(f'ERROR: {len(db)} signatures found. Need exactly one.',
              file=sys.stderr)
        sys.exit(-1)

    items = mh.hashes
    chooser = list(items.keys())
    tracking = {}
    tracking_mh = mh.copy_and_clear()

    for k in chooser:
        tracking[k] = items[k]

    xx = []
    n = 0
    while chooser:
        for i in range(args.sampling_interval):
            h = random.choice(chooser)
            val = tracking.get(h)
            if val is not None:
                if val <= 0:
                    del tracking[h]
                    continue

                n += 1
                tracking[h] -= 1
                tracking_mh.add_hash(h)

        hh = tracking_mh.hashes
        gt2 = 0
        for k, v in hh.items():
            if v > 1:
                gt2 += 1
        
        sat = gt2 / len(tracking_mh)

        chooser = list(tracking)
        xx.append((n, len(tracking_mh), tracking_mh.sum_abundances, sat))
        print(xx[-1])

    with open(args.output_csv, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['n', 'intersect_bp', 'weighted_bp', 'saturation'])

        scaled = args.scaled
        for (n, isect, wt, sat) in xx:
            w.writerow([n, isect*scaled, wt*scaled, sat])


if __name__ == '__main__':
    sys.exit(main())
