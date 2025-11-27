import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the scores from MrOS, SOL, Framingham Heart Study, and Rotterdam Study
scores_server1 = pd.read_csv('combined_pca_scores.csv', index_col = 0)
#Load the scores from the Doetinchem Cohort Study
scores_dcs = pd.read_csv('pca_scores_DCS_common_genera.csv', index_col = 0)
#Load the scores from LifeLines
scores_ll = pd.read_csv('pca_scores_DAGIII_common_genera.csv', index_col = 0)
scores_ll['Cohort'] = 'LL'

#Combine scores into one dataframe
scores_combined = scores_server1.merge(scores_dcs, how = 'outer').merge(scores_ll, how = 'outer')

#Reshuffle to avoid order of combining cohort determines plot
scores = scores_combined.sample(frac=1)

#Load explained variance per PCA
explained_variance = pd.read_csv('explained_variance.csv', index_col = 0)

# Plot PCA
sns.set(style="whitegrid")
plt.figure(figsize=(10, 8))
ax = sns.scatterplot(data=scores, x="PC1", y="PC2", hue="Cohort", palette="Set2", s=75, alpha=0.5)

# Sort legend alphabetically s
handles, labels = ax.get_legend_handles_labels()
# First entry is the legend title; keep it separate
handles = handles[1:]
labels = labels[1:]
labels_handles_sorted = sorted(zip(labels, handles), key=lambda x: x[0])
sorted_labels, sorted_handles = zip(*labels_handles_sorted)
plt.legend(sorted_handles, sorted_labels, title="Cohort", loc="best", fontsize=10)
plt.title("PCA of CLR-Transformed Data for Top 7 Genera (Grouped by Cohort)", fontsize=14)
plt.xlabel(f"PC1 ({explained_variance['0'][0]*100:.1f}%)")
plt.ylabel(f"PC2 ({explained_variance['0'][1]*100:.1f}%)")
plt.tight_layout()
plt.savefig("combined_pca_plot.png")
plt.show()
