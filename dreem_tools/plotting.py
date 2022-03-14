import matplotlib.pyplot as plt
import seaborn as sns

def colors_for_sequence(seq : str):
    colors = []
    for e in seq:
        if e == 'A':
            colors.append('red')
        elif e == 'C':
            colors.append('blue')
        elif e == 'G':
            colors.append('orange')
        else:
            colors.append('green')
    return colors


def plot_pop_avg(seq, ss, reactivities, ax=None):
    colors = colors_for_sequence(seq)
    x = list(range(len(seq)))
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(20, 4))
    ax.bar(range(0, len(reactivities)), reactivities, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n{nt}" for s, nt in zip(seq, ss)])
    return ax

