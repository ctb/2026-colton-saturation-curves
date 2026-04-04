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

    #
    # load the signature
    #

    db = sourmash.load_file_as_index(args.sketch)
    db = db.select(ksize=args.ksize, scaled=args.scaled)
    if len(db) != 1:
        print(f'ERROR: {len(db)} signatures found. Need exactly one.',
              file=sys.stderr)
        sys.exit(-1)

    ss = list(db.signatures())[0]
    mh = ss.minhash.downsample(scaled=args.scaled)
    print(f"Loaded signature '{ss.name}' - len {len(mh)} / sum {mh.sum_abundances}")
    del ss

    if not mh.track_abundance:
        print(f'ERROR: {len(db)} signatures found. Need exactly one.',
              file=sys.stderr)
        sys.exit(-1)

    #
    # extract the hashes with their abundances
    #

    items = mh.hashes
    chooser = list(items.keys())
    tracking = {}
    sampling_mh = mh.copy_and_clear()

    # copy things over to a modifiable dictionary { hash: abundance }
    for k in chooser:
        tracking[k] = items[k]

    #
    # iterate at sampling interval; at each iteration,
    # * randomly choose hashes
    # * get abundance from 'tracking' and subtract one
    # * if abundance drops below 0, remove;
    # 

    samples = []
    n = 0
    while chooser:
        for i in range(args.sampling_interval):
            h = random.choice(chooser)
            val = tracking.get(h)
            if val is None:     # no longer present
                continue

            assert val > 0

            # decrement; remove if 0
            if val == 1:
                del tracking[h]
            else:
                tracking[h] -= 1

            # track sampled abund.
            n += 1
            sampling_mh.add_hash(h)

        # at each sampling interval, track how many abund 1, 2, 3, 4, 5.
        hh = sampling_mh.hashes
        gt_d = {}
        for i in range(1, 6):
            gt_d[i] = 0

        gt1 = 0
        
        for h, abund in hh.items():
            if abund > 1:
                gt1 += 1
            for i in range(1, 6):
                if abund >= i:
                    gt_d[i] += 1

        assert len(gt_d) <= 5   # 1-5, not == 0 or > 5

        yy = []
        for abund in range(1, 6):
            frac = gt_d[abund] / len(sampling_mh)
            yy.append(frac)

        assert len(yy) == 5

        sat2 = gt1 / len(sampling_mh)
        
        samples.append([n, len(sampling_mh), sampling_mh.sum_abundances, sat2]\
                       + yy)
        print(samples[-1])
        assert len(samples[-1]) == 9

        # reinitialize key choices from tracking list
        chooser = list(tracking)

    # done!

    with open(args.output_csv, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['n', 'intersect_bp', 'weighted_bp', 'sat',
                    'f1', 'f2', 'f3', 'f4', 'f5'])

        scaled = args.scaled
        for (n, isect, wt, sat, *rest) in samples:
            w.writerow([n, isect*scaled, wt*scaled, sat] + rest)


if __name__ == '__main__':
    sys.exit(main())
