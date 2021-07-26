import matplotlib.pyplot as plt
import numpy as np


def volcano(data, subsets, ax=None, interactive=False):
    """
    make a volcano plot

    :param data: dataframe containing comparison data
    :param subsets: list of subsets of points to have different formatting, as tuples (set, plotting kwargs, lookup)
    :param ax: axes to plot on
    """
    layers = [
    ]  # list of layers of points (set, plotting kwargs) with distinct points
    prev = set()  # set of points already plotted on previous layers
    for s in subsets[::-1]:  # iterate from last subset backwards
        if s[2] is None:
            # subtract points already plotted on previous layers
            layers.append((list(s[0]-prev), s[1]))
        else:
            layers.append((
                list(data[data[s[2]].apply(lambda x: type(x) == str and (
                    x in s[0]-prev or any([a in s[0]-prev for a in x.split(";")])))].index),
                s[1]
            ))  # lookup points to get index, then subtract points already plotted on previous layers
        prev |= s[0]  # add current layer's points to prev
    layers = layers[::-1]  # reverse layers

    if ax is None:
        ax = plt.gca()
    # iterate from last layer backwards so that the last subset is plotted on top
    paths = []
    for s in layers:
        path = ax.scatter(data.loc[s[0]]["log FC"], -
                          np.log10(data.loc[s[0]]["p adjusted"]), **s[1])
        paths.append(path)

    if interactive:
        # handle interactive annotations
        fig = ax.get_figure()
        # make dummy annotation
        annot = ax.annotate("", xy=(0, 0), xytext=(10, 20), textcoords="offset points", ha="center",
                            bbox=dict(boxstyle="round", fc="w", alpha=0.4),
                            arrowprops=dict(arrowstyle="fancy", connectionstyle="arc3,rad=0.2", fc="k", ec="k"))
        annot.set_visible(False)

        def update_annot(i, ind):
            # show annotation for ind-th point of i-th layer
            idx = layers[i][0][ind]  # get index of point in data
            annot.xy = (data.loc[idx]["log FC"], -
                        np.log10(data.loc[idx]["p adjusted"]))
            text = "{}\n{:.2f}, {:.2f}".format(
                data.loc[idx]["gene"],
                data.loc[idx]["log FC"],
                -np.log10(data.loc[idx]["p adjusted"]))
            annot.set_text(text)

        def hover(event):
            if event.inaxes == ax:
                # get mouse coordinates normallized to axis range
                xrange, yrange = ax.get_xlim(
                )[1] - ax.get_xlim()[0], ax.get_ylim()[1] - ax.get_ylim()[0]
                x, y = event.xdata / xrange, event.ydata / yrange

                clear = True  # clear annotation if not hovering over a point
                for i, path in list(enumerate(paths))[::-1]:
                    cont, ind = path.contains(event)
                    if cont:
                        clear = False
                        # calculate distances to each point mouse is hovering over
                        distances = [
                            np.linalg.norm(
                                np.array([x, y]) - np.array(
                                    [data.loc[layers[i][0][z]]["log FC"] / xrange,
                                     -np.log10(data.loc[layers[i][0][z]]["p adjusted"]) / yrange]
                                )
                            ) for z in ind["ind"]
                        ]
                        # get index of closest point
                        index = ind["ind"][np.argmin(distances)]
                        # update annotation text for the closest point
                        update_annot(i, index)
                        annot.set_visible(True)
                        fig.canvas.draw_idle()
                        break
                if clear:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)
