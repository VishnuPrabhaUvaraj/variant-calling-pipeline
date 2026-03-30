import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip

print("Starting variant analysis...")

# ── 1. Load and parse VCF ──────────────────────────
def parse_vcf(vcf_path):
    rows = []
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            chrom, pos, _, ref, alt, qual, filt = cols[:7]
            vtype = 'SNP' if len(ref)==1 and len(alt)==1 else 'INDEL'
            rows.append({
                'CHROM': chrom,
                'POS': int(pos),
                'REF': ref,
                'ALT': alt,
                'QUAL': float(qual) if qual != '.' else 0,
                'FILTER': filt,
                'TYPE': vtype
            })
    return pd.DataFrame(rows)

df = parse_vcf('results/variants/NA12878.final.vcf.gz')
df_pass = df[df['FILTER'] == 'PASS']

print(f"Total variants:  {len(df)}")
print(f"PASS variants:   {len(df_pass)}")
print(f"SNPs:            {(df_pass.TYPE=='SNP').sum()}")
print(f"INDELs:          {(df_pass.TYPE=='INDEL').sum()}")

# ── 2. Figure 1+2+3 in one image ───────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Pie chart
type_counts = df_pass['TYPE'].value_counts()
axes[0].pie(
    type_counts,
    labels=type_counts.index,
    autopct='%1.1f%%',
    colors=['#378ADD', '#1D9E75'],
    startangle=90,
    wedgeprops={'linewidth': 1, 'edgecolor': 'white'}
)
axes[0].set_title('Variant types (PASS)', fontsize=13)

# Quality score histogram
axes[1].hist(df_pass['QUAL'], bins=50, color='#7F77DD',
             edgecolor='white', linewidth=0.5)
axes[1].set_xlabel('Quality score (QUAL)')
axes[1].set_ylabel('Count')
axes[1].set_title('Variant quality distribution', fontsize=13)
axes[1].axvline(30, color='#D85A30', linewidth=1.5,
                linestyle='--', label='Q30')
axes[1].legend()

# Ti/Tv ratio
ti_pairs = {('A','G'),('G','A'),('C','T'),('T','C')}
snps = df_pass[df_pass['TYPE']=='SNP'].copy()
snps['TiTv'] = snps.apply(
    lambda r: 'Ti' if (r['REF'], r['ALT']) in ti_pairs else 'Tv',
    axis=1
)
ti = (snps['TiTv']=='Ti').sum()
tv = (snps['TiTv']=='Tv').sum()
ratio = ti / tv if tv > 0 else 0
titv_counts = pd.Series({'Transitions (Ti)': ti, 'Transversions (Tv)': tv})
axes[2].bar(titv_counts.index, titv_counts.values,
            color=['#378ADD', '#D85A30'], width=0.5)
axes[2].set_title(f'Ti/Tv ratio: {ratio:.2f}', fontsize=13)
axes[2].set_ylabel('Count')
for i, v in enumerate(titv_counts.values):
    axes[2].text(i, v + 10, str(v), ha='center', fontsize=11)

plt.tight_layout()
plt.savefig('figures/variant_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: figures/variant_summary.png")
print(f"Ti/Tv ratio: {ratio:.3f}  (expected ~2.0-2.5 for WES)")

# ── 3. Figure 4: Coverage depth plot ───────────────
print("Loading coverage data (this may take a moment)...")
cov = pd.read_csv(
    'results/alignment/NA12878.coverage.txt',
    sep='\t', header=None,
    names=['chrom', 'pos', 'depth']
)
cov_sample = cov.sample(min(10000, len(cov)), random_state=42)
cov_sample = cov_sample.sort_values('pos')

fig2, ax = plt.subplots(figsize=(14, 4))
ax.fill_between(cov_sample['pos'], cov_sample['depth'],
                alpha=0.4, color='#1D9E75')
ax.plot(cov_sample['pos'], cov_sample['depth'],
        color='#0F6E56', linewidth=0.5)
ax.axhline(
    cov['depth'].mean(),
    color='#D85A30',
    linewidth=1.5,
    linestyle='--',
    label=f"Mean: {cov['depth'].mean():.1f}x"
)
ax.set_xlabel('Position on chr20')
ax.set_ylabel('Depth of coverage')
ax.set_title('Coverage depth across chr20', fontsize=13)
ax.legend()
plt.tight_layout()
plt.savefig('figures/coverage_plot.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: figures/coverage_plot.png")

# ── 4. Effect type summary ─────────────────────────
print("\n=== Variant Effect Summary ===")
annot_file = 'results/annotation/NA12878.variant_table.tsv'
effects = []
with open(annot_file) as f:
    for line in f:
        for effect in ['missense_variant', 'synonymous_variant',
                       'intron_variant', 'upstream_gene_variant',
                       'downstream_gene_variant']:
            if effect in line:
                effects.append(effect)
                break

effect_counts = pd.Series(effects).value_counts()
print(effect_counts.to_string())

fig3, ax = plt.subplots(figsize=(10, 5))
effect_counts.plot(kind='bar', color='#7F77DD', ax=ax)
ax.set_xlabel('Effect type')
ax.set_ylabel('Count')
ax.set_title('Variant effects distribution', fontsize=13)
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.savefig('figures/effect_distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: figures/effect_distribution.png")

print("\nAll figures saved to figures/")
print("Phase 7 complete!")
