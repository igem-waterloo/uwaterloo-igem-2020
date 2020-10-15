import os
import re
import numpy as np
import matplotlib.pyplot as plt


def generate_heatmap(mode="none", log_scale=False, eps=1e-2):

    """
    Displays a heatmap given a mode indicating how to transform the data


    :param mode: How the data should be processed. Options are none, relative, scaled, truncated
    none - Do not transform data at all
    relative - Subtract unmutated scores from all results
    scaled - Scales all columns in [0, 1]
    truncated - Same as none except overall energy scores are removed
    :param log_scale: Take logarithm of all values
    :param eps: Added to all values during log scaling to prevent dividion by zero. Is ignored otherwise
    :return: None
    """

    # This code assumes that every file in the directory given is a pdb file. Create a directory and fill it with
    # the files in need of processing, point 'base' at the directory, and run the script.

    mode = mode.lower()

    assert mode == 'none' or mode == 'relative' or mode == 'scaled' or mode == 'truncated'

    base = 'Rosetta\\mutated_pdb'  # this needs to be the directory containing all the data and nothing else
    direct = os.listdir(base)
    energy = []
    resnames = []

    for file in direct:
        if file.endswith('pdb'):
            path = os.path.join(base, file)
            pdb = open(path)
            beng = []
            bres = []
            intable = False
            for line in pdb:

                if "#BEGIN_POSE_ENERGIES_TABLE" in line:
                    intable = True
                    continue

                if "#END_POSE_ENERGIES_TABLE" in line:
                    intable = False

                if intable:
                    line = re.sub(r"\s+", " ", line)
                    bres.append(line.split(" ")[:1])
                    fval = line.split(" ")[-2:][0]
                    try:
                        fval = float(fval)
                        beng.append(fval)
                    except Exception:
                        pass

            energy.append(beng)

    bres = bres[(len(bres) - len(energy[0])):]

    energy = np.array(energy)

    if mode == 'truncated':
        energy = energy.transpose()
        energy = energy[1:]
        energy = energy.transpose()
        bres = bres[1:]
    if mode == 'relative':
        energy = energy - (energy[-1:])
    if mode == 'scaled':
        energy = (energy - energy.min(0)) / (energy.max(0) - energy.min(0))

    if log_scale:
        energy = np.log(energy - energy.min() + eps)

    figure, axis = plt.subplots()
    im = axis.imshow(energy)

    xlabels = []
    for v in bres:
        xlabels.append(v[0])

    axis.set_xticks(np.arange(len(bres)))
    axis.set_xticklabels(xlabels)

    ylabels = []
    for v in direct:
        if v.endswith('pdb'):
            ylabels.append(v.replace('HETATM_relaxed', ''))

    axis.set_yticks(np.arange(len(ylabels)))
    axis.set_yticklabels(ylabels)

    axis.set_ylabel("Proteins")
    axis.set_xlabel("Residues")

    plt.title('Per Residue Energy Scores for Mutant and Original Proteins')

    cbar = plt.colorbar(im)
    desc = ""
    if mode == "scaled":
        desc = "Scaled "
    elif mode == "relative":
        desc = "Relative "
    cbar.ax.set_ylabel("{}Energy Score {}(Rosetta Energy Units)".format(desc, 'log' if log_scale else ''))

    plt.setp(axis.get_xticklabels(), rotation=90, ha="right", va="center", rotation_mode="anchor")
    plt.tight_layout()

    plt.show()


if __name__ == '__main__':

    mode = input("Input heatmap mode (none, scaled, relative, truncated):")
    generate_heatmap(mode)

