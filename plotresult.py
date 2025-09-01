import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
import arviz as az

import corner
import pickle

import seaborn as sns
import pandas as pd
from collections import OrderedDict
import elfi

rcParams["font.size"] = 22
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

# Load the result_BOLFI object from the file
with open("sample_nsample100000_result_bolfi_fail_detecton_multi_euclidean_data_n1000_Nov14.pkl", "rb") as f:
    result_BOLFI = pickle.load(f)

###########################
#data = pd.DataFrame(result_BOLFI.samples)

# Create a corner plot using sns pairplot
#sns.pairplot(data, diag_kind='kde', markers='o', plot_kws={'alpha': 0.5})

#plt.suptitle('BOLFI sample Corner Plot Example', y=1.02, fontsize=16)
# Save the current figure
#plt.savefig("sample_nsample100000_result_bolfi_fail_detecton_multi_euclidean_data_n1000_Nov14_cornerplot.jpg", format="jpg", bbox_inches="tight",dpi=150)



###########################
#using corner package
# Set up the parameters of the problem.
ndim, nsamples = 6, 200000



data= dict(result_BOLFI.samples)
# Create a [6, N] array from the data
# We will stack the arrays along the first axis
data_array = np.array([data['tauE'],
                       data['tauI'],
                       data['prna'],
                       data['p'],
                       data['gamma'],
                       data['beta']]).T

# Plot it.
figure = corner.corner(data_array, labels=[r"$\tau_{E}$", r"$\tau_{I}$", r"$\log_{10}{Prna}$", r"$\log_{10}{p}$",
                                           r"$\log_{10}{\gamma}$", r"$\log_{10}{\beta}$"],
                       quantiles=[0.025, 0.5, 0.975],
                       show_titles=True, title_kwargs={"fontsize": 16})
figure.gca().annotate("BOLFI result corner plot",
                      xy=(1.0, 1.0), xycoords="figure fraction",
                      xytext=(-20, -10), textcoords="offset points",
                      ha="right", va="top")
figure.savefig("BOLFI_cornerReorder.png", dpi=300)

