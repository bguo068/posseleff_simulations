import matplotlib.pyplot as plt
import numpy as np


def area(a):
    """
    calculate area under curve y = x^a between[0, 1]
    """
    return 1 / (a + 1)


def x_a(x):
    return np.power(x, a)


def y_orig(x: np.ndarray):
    """
    if (!sim_relatedness)
            return T;
    if (child.index % 2 == 0)
            return T;

    prob = parent1.relatedness(parent2);
    if ((prob > 1.0 / 32) & (prob < 1.0 / 4))
        return T;
    else
        return sample(c(T, F), size=1, replace=T, weights=c(1, 999));
    """
    sel = (x > 1.0 / 32) & (x < 1.0 / 4)
    y = np.ones(x.shape)
    y[sel] = 0.5 + 0.5
    y[~sel] = 0.5 + 0.5 * 0.001
    return y


def area_orig():
    L = (1 / 4) - (1 / 32)
    area = L * (0.5 + 0.5) + (1 - L) * (0.5 + 0.5 * 0.001)
    return area


fig, ax = plt.subplots()
for a in [0, 1 / 10, 1 / 3, 1, 3, 10, "current"]:
    x = np.arange(0, 1, 0.01)
    if a == "current":
        y = y_orig(x) / area_orig()
        ax.plot(x, y, "k", lw=3, label="current")
    else:
        y = x_a(x) / area(a)
        if a != 0:
            ax.plot(x, y, label=f"power={a:.2f}")
        else:
            ax.plot(x, y, "r--", lw=2, label=f"flat")

    ax.set_xlabel("relatedness")
    ax.set_ylabel("acceptance probability")
ax.legend()

fig.savefig("relatedness_transform.pdf")

