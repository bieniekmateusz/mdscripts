#!/usr/bin/python
"""
Based modelled on Enrico's script "calculation-s2-dipole-backbone-aa.py". The calculations have been partly checked against Enrico's results which in turn have been tested earlier.

Calculation of s2 of the N-H of the backbone using the equation 12 of the following article
http://www.sciencedirect.com/science/article/pii/S0022283684700908

Mateusz Bieniek & Enrico Spiga
"""

import sys
import os
import numpy as np
import MDAnalysis
import MDAnalysis.analysis.distances
import argparse


def s2(universe, frames_no=-1):
    # extract N-H bonds
    hn = universe.select_atoms('name HN and not resname PRO')
    n = universe.select_atoms('name N and not resname PRO and not bynum 1')    # ignore the first N

    # for each residue there are two atoms, which are just after the other
    # e.g. atom 20 corresponds to 21, atoms 41 corresponds to 42, etc
    assert all([left.index == right.index + 1 for left, right in zip(hn, n)])

    # corresponding resids of the selected atoms
    resids = hn.resids

    # create zero matrices 3x3 for each residue
    # will hold the sums for each residue
    sums_3X3 = [np.zeros([3, 3]) for _ in resids]

    # all frames ?
    if frames_no == -1: frames_no = len(u.trajectory)

    for _ in u.trajectory[:frames_no]:
        # the [2] extracts the distances (the first two are residue number of the atoms)
        dnorms = MDAnalysis.analysis.distances.dist(hn, n)[2]

        # N-H dists in each dimension
        dists = n.atoms.positions - hn.atoms.positions

        # normalise (divide the dists in each dimension by their unit vector)
        for i in range(len(dnorms)):
            normalised_dist = dists[i] / dnorms[i]
            # convert the numpy matrix
            m_dists = np.matrix(normalised_dist)
            # multiply each by each
            each_by_each = np.transpose(m_dists) * m_dists
            sums_3X3[i] += each_by_each


    # calculate the s2 order parameter
    s2s = []
    for summed_matrix in sums_3X3:
        # average
        summed_matrix /= frames_no
        summed_matrix **= 2
        np_sum = np.sum(summed_matrix)

        s2 = 0.5 * (3 * np_sum - 1)
        s2s.append(s2)

    # all values have to be between 0 and 1
    assert all([0 <= s2 <= 1 for s2 in s2s])

    # add the residue indexes to each s2
    s2s_res = zip(resids, s2s)

    return s2s_res


if __name__ == "__main__":
    # argument parser
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-s", help="path to topology file", metavar="topology", required=True)
    argparser.add_argument("-f", help="path to trajectory file", metavar="trajectory", required=True)
    argparser.add_argument("-o", help="path to the output file with s2 order parameters", metavar="output", required=True)
    argparser.add_argument("-e", help="the last frame for the analaysis", metavar="last_frame", type=int)

    # parse arguments
    args = sys.argv[1:]
    args = argparser.parse_args(args)
    topo_s = args.s
    traj_f = args.f
    output_s2_file = args.o
    last_frame = args.e

    # was last frame given?
    if last_frame is None:  last_frame = -1

    # files?
    assert os.path.isfile(topo_s)
    assert os.path.isfile(traj_f)

    # do not let overwrite
    if os.path.isfile(output_s2_file):
        print 'I will not overwrite this file! Leave "%s" in peace. ' % output_s2_file
        sys.exit(1)

    # read the trajectory using MDAnalysis
    u = MDAnalysis.Universe(topo_s, traj_f)
    print "Number of all frames in the trajctory: ", len(u.trajectory)

    # compute s2 order parameters
    s2s_data = s2(u, frames_no=last_frame)

    # output to file
    output = os.linesep.join(["%d %f" % (resid, nhS2s) for resid, nhS2s in s2s_data])
    open(output_s2_file, 'w').write(output)

    print
    print "The s2 order parameters calculated successfully. Hopefully. "