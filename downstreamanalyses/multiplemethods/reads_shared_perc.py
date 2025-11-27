import pandas as pd

# Define file paths
cohort = "FHS"
file_16s = f"{cohort}_16S_Prev_Abundance_per_genus.csv"
file_wgs = f"{cohort}_WGS_Prev_Abundance_per_genus.csv"

# Load data
df_16s = pd.read_csv(file_16s)
df_wgs = pd.read_csv(file_wgs)

# Drop rows without genus info (optional)
df_16s = df_16s[df_16s['Genus'].notna()]
df_wgs = df_wgs[df_wgs['Genus'].notna()]

# Unique genera
genera_16s = set(df_16s['Genus'])
genera_wgs = set(df_wgs['Genus'])

# Shared genera
shared_genera = genera_16s & genera_wgs

# Counts
n_16s = len(genera_16s)
n_wgs = len(genera_wgs)
n_shared = len(shared_genera)

# Percent shared
pct_16s = n_shared / n_16s * 100
pct_wgs = n_shared / n_wgs * 100

# % of reads (mean abundance) covered by shared genera
reads_16s = df_16s[df_16s['Genus'].isin(shared_genera)]['Mean_Abundance'].sum() * 100  # relative freq to %
reads_wgs = df_wgs[df_wgs['Genus'].isin(shared_genera)]['Mean_Abundance'].sum() * 100

# Output
print(f"In {cohort}, {n_shared} genera were shared between sequencing methods,")
print(f"representing {pct_16s:.1f}% of {n_16s} 16S genera and {pct_wgs:.1f}% of {n_wgs} shotgun genera.")
print(f"These shared genera explained {reads_16s:.1f}% of 16S reads and {reads_wgs:.1f}% of shotgun reads.")
