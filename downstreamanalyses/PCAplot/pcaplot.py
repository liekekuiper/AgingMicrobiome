import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

scores_server1 = pd.read_csv('combined_pca_scores.csv', index_col = 0)
scores_dcs = pd.read_csv('pca_scores_DCS_common_genera.csv', index_col = 0)
scores = scores_server1.merge(scores_dcs, how = 'outer')
explained_variance = pd.read_csv('explained_variance.csv', index_col = 0)

# Plot PCA
sns.set(style="whitegrid")
plt.figure(figsize=(10, 8))
sns.scatterplot(data=scores, x="PC1", y="PC2", hue="Cohort", palette="Set2", s=100, alpha=0.5)
plt.title("PCA of CLR-Transformed Data for Top 4 Genera (Grouped by Cohort)", fontsize=14)
plt.xlabel(f"PC1 ({explained_variance['0'][0]*100:.1f}%)")
plt.ylabel(f"PC2 ({explained_variance['0'][1]*100:.1f}%)")
plt.legend(title="Cohort", loc="best", fontsize=10)
plt.tight_layout()
plt.savefig("combined_pca_plot.png")
plt.show()
