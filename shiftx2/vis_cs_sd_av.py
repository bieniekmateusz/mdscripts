#!/usr/bin/python

"""
    Quick and dirty script. Dirty dirty dirty. 
    Visualise the output of the 'shiftx2_cs_sd_av' script.
"""

import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# for conversion
threeToOne = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
oneToThree = dict(zip(threeToOne.values(), threeToOne.keys()))

def read_cs_from_file(filename):
    lines = open(filename).read().split(os.linesep)
    split = [line.split() for line in lines if line.strip() is not '']
    data = [(int(resid), res_name_1_letter, float(cs_av), float(cs_sd)) for resid, res_name_1_letter, cs_av, cs_sd in split]
    # well, apparently one too many ...
    return data[:-1]

# read the charm data
data = read_cs_from_file('/home/dresio/phd/chemical_shift_comparison/copiedCharmAvSd/cs_sum_dBr1c36_every10th')
res_num = zip(*data)[0]
res_aaa = [oneToThree[res_letter] for res_letter in zip(*data)[1]]
cs_av = zip(*data)[2]
cs_sd = zip(*data)[3]
res_num_from_1 = list(range(1, len(cs_av) + 1))

# read the cs extracted from nmr
nmr_cs_filename = '/home/dresio/phd/chemical_shift_comparison/copiedCharmAvSd/extracted4283chemicalShiftsN15_firstNs'
lines = open(nmr_cs_filename).read().split(os.linesep)
split = [line.split() for line in lines]
nmr_cs_data = [(int(resid), three_letter_aa, float(cs)) for resid, three_letter_aa, _, _, _, cs in split]
nmr_res_nums = zip(*nmr_cs_data)[0]
nmr_res_aaa = zip(*nmr_cs_data)[1]
nmr_cs = zip(*nmr_cs_data)[2]

# the names correspond to each other?
assert all([x == y for x,y in zip(res_aaa, nmr_res_aaa)])
assert len(res_aaa) == len(nmr_res_aaa)

for x, nmr_cs_point, cs_av_point in zip(res_num_from_1, nmr_cs, cs_av):
    plt.plot((x, x), (nmr_cs_point, cs_av_point), 'k--', linewidth=0.5, color='grey')

# plot charm
plt.errorbar(res_num_from_1, cs_av, yerr=cs_sd, fmt='o', label='Charmm 36')
# plot nmr
plt.plot(res_num_from_1, nmr_cs, 'bs', color='r', label='NMR')

# plot x label & number residues from 1
enumerated_ticks = [str(index) + ' ' + aaa for index, aaa in enumerate(res_aaa, start=1)]
plt.xticks(res_num_from_1, enumerated_ticks, rotation='vertical', family='monospace')
# do not let the points overlap with the y axis
plt.xlim([res_num_from_1[0] - 1, res_num_from_1[-1] +1])

# add the best match
# lines = open('/home/dresio/phd/chemical_shift_comparison/copiedCharmAvSd/converted9418.pdb.cs-bestmatch').read().split('\n')
# #lines = open('/home/dresio/phd/chemical_shift_comparison/copiedCharmAvSd/best_nonsq.pdb.cs').read().split('\n')
# lines = [line.split(',') for line in lines][:-1]
# best_match_data = [(int(resid), a, atom, float(cs)) for resid, a, atom, cs in lines]
# best_match_data = filter(lambda x: x[2] == 'N', best_match_data)
# # remove ala
# best_match_data = best_match_data[:-1]
# best_match_cs = zip(*best_match_data)[-1]
# plt.plot(res_num_from_1, best_match_cs, 'bs', color='green')


# draw the secondary structure
"""
The following is the secondary structure as downloaded from the web page:
QLTDLSFVDITDSSIGLRWTLNSSTIIGYRITVVAAGEGIIFEDFVDSSVGYYTVTGLEGIDYDISVITLINGGESATTLTQQT
  TT EEE   SS EEEE     SS  EEEEEEEEETS  EEEEE  SS SEEEE    TSEEEEEEEEEETTEE   EEEEE

Where E:beta strand, T:turn, S: bend

Note 1: that the extra residues on both ends have been deleted.
Note 2: proline residues and their corresponding classificaiton has been removed
"""
# ntoe that the prolines have been deleted ...
secondary_strs = '  TT EEE   SS EEEE     SS  EEEEEEEEETS  EEEEE  SS SEEEE    TSEEEEEEEEEETTEE   EEEEE '
assert len(secondary_strs) == len(res_aaa)
currentAxis = plt.gca()
for i, res in enumerate(secondary_strs):
    height = 2
    if res == 'E':
        # beta strand
        rect = plt.Rectangle((i,100), 1, height, color='yellow')
        currentAxis.add_patch(rect)
    elif res == 'T':
        # turn
        rect = plt.Rectangle((i,100), 1, height, color='purple', alpha=1)
        currentAxis.add_patch(rect)
    elif res == 'S':
        # bend
        rect = plt.Rectangle((i,100), 1, height, color='blue', alpha=1)
        currentAxis.add_patch(rect)

plt.legend(loc='upper left')
plt.ylabel('Nitrogen-15')

# add the rmsf
lines = open('/home/dresio/phd/chemical_shift_comparison/copiedCharmAvSd/dBr1c36/rmsf-renum-ns.xvg-justdata').read().split(os.linesep)
lines = filter(None, [line.split() for line in lines])
rmsf = [(res_num, rmsf) for res_num, rmsf in lines]
rmsf_val = [float(val) for val in list(zip(*rmsf)[1])]

assert 'PRO' not in res_aaa

rmsf_res = 'GLN LEU THR ASP LEU SER PHE VAL ASP ILE THR ASP SER SER ILE GLY LEU ARG TRP THR PRO LEU ASN SER SER THR ILE ILE GLY TYR ARG ILE THR VAL VAL ALA ALA GLY GLU GLY ILE PRO ILE PHE GLU ASP PHE VAL ASP SER SER VAL GLY TYR TYR THR VAL THR GLY LEU GLU PRO GLY ILE ASP TYR ASP ILE SER VAL ILE THR LEU ILE ASN GLY GLY GLU SER ALA PRO THR THR LEU THR GLN GLN THR ALA'.split()
assert len(rmsf_res) == len(rmsf)
# find the prolines (PRO)
prolines = [index for index, aaa in enumerate(rmsf_res) if aaa == 'PRO']
# remove prolines and their corresponding values
for pro_index in prolines[::-1]:
    del rmsf_res[pro_index]
    del rmsf_val[pro_index]

# remove the last ALA residue from rmsf
del rmsf_res[-1]
del rmsf_val[-1]

assert len(rmsf_val) == len(res_aaa)
assert rmsf_res == res_aaa

t = plt.twinx()
t.plot(res_num_from_1, rmsf_val, label='rmsf', alpha=0.8, color='black')
plt.ylabel(r'RMSF ($\AA$)')

plt.title('Comparison of chemical shifts from simulation and NMR')
plt.xlabel('Residue')
plt.xlim([res_num_from_1[0] - 1, res_num_from_1[-1] +1])

plt.legend(loc='upper right')

# find correlation between the simulation chemical shift standard deviation and rmsf
import numpy as np
# print np.cov(cs_sd, rmsf_val)
print 'correlation sd & val', np.corrcoef(cs_sd, rmsf_val)
print 'corrr av & nmr', np.corrcoef(cs_av, nmr_cs)
# plt.show()
