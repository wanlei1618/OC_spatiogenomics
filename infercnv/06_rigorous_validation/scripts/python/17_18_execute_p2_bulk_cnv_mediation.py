import csv
import math
import os
import random
import zipfile
from collections import defaultdict
from datetime import datetime
from xml.sax.saxutils import escape


PROJECT_ROOT = r"D:\OC_spatiogenomics\infercnv"
OUT_ROOT = os.path.join(PROJECT_ROOT, "06_rigorous_validation")
SRC_BULK = os.path.join(PROJECT_ROOT, "03_spp1_cd44_itgb1_hypothesis_validation", "SPP1_ITGB1_CD44_hypothesis_validation_complete", "tables")
SRC_CNV = os.path.join(PROJECT_ROOT, "CNV_expression_joint_analysis_integrated_oc", "tables")

DIRS = {
    "bulk": os.path.join(OUT_ROOT, "11_bulk_negative_validation"),
    "evidence": os.path.join(OUT_ROOT, "12_integrated_evidence_scoring"),
    "tables": os.path.join(OUT_ROOT, "tables"),
    "reports": os.path.join(OUT_ROOT, "reports"),
    "logs": os.path.join(OUT_ROOT, "logs"),
}
for path in DIRS.values():
    os.makedirs(path, exist_ok=True)

RANDOM_SEED = 20260711
random.seed(RANDOM_SEED)
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_PATH = os.path.join(DIRS["logs"], f"17_18_execute_p2_bulk_cnv_mediation_{STAMP}.log")


def read_csv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fnum(value):
    try:
        if value in ("", None):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def clean_pairs(x, y):
    return [(a, b) for a, b in zip(x, y) if not math.isnan(a) and not math.isnan(b)]


def mean(values):
    values = [v for v in values if not math.isnan(v)]
    return sum(values) / len(values) if values else math.nan


def median(values):
    values = sorted(v for v in values if not math.isnan(v))
    if not values:
        return math.nan
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def rankdata(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = rank
        i = j + 1
    return ranks


def pearson(x, y):
    pairs = clean_pairs(x, y)
    n = len(pairs)
    if n < 3:
        return math.nan, n
    xs, ys = zip(*pairs)
    mx, my = mean(xs), mean(ys)
    sx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    sy = math.sqrt(sum((b - my) ** 2 for b in ys))
    if sx == 0 or sy == 0:
        return math.nan, n
    r = sum((a - mx) * (b - my) for a, b in pairs) / (sx * sy)
    return r, n


def normal_p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def spearman(x, y):
    pairs = clean_pairs(x, y)
    n = len(pairs)
    if n < 3:
        return math.nan, math.nan, n
    xs, ys = zip(*pairs)
    r, _ = pearson(rankdata(list(xs)), rankdata(list(ys)))
    if math.isnan(r) or abs(r) >= 1 or n <= 3:
        return r, math.nan, n
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(max(n - 3, 1))
    return r, normal_p_from_z(z), n


def bh_adjust(rows, p_col="p_value", out_col="p_adj_BH"):
    valid = sorted([(i, fnum(row.get(p_col))) for i, row in enumerate(rows) if not math.isnan(fnum(row.get(p_col)))], key=lambda x: x[1])
    m = len(valid)
    adjusted = [math.nan] * len(rows)
    running = 1.0
    for rank, (i, p) in reversed(list(enumerate(valid, start=1))):
        running = min(running, p * m / rank)
        adjusted[i] = running
    for i, row in enumerate(rows):
        row[out_col] = adjusted[i]


def random_effects_meta(effects, variances):
    data = [(e, v) for e, v in zip(effects, variances) if not math.isnan(e) and not math.isnan(v) and v > 0]
    k = len(data)
    if k == 0:
        return {"k": 0}
    eff = [e for e, _ in data]
    var = [v for _, v in data]
    w = [1 / v for v in var]
    fixed = sum(wi * ei for wi, ei in zip(w, eff)) / sum(w)
    q = sum(wi * (ei - fixed) ** 2 for wi, ei in zip(w, eff))
    c = sum(w) - sum(wi ** 2 for wi in w) / sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = [1 / (v + tau2) for v in var]
    pooled = sum(wi * ei for wi, ei in zip(wr, eff)) / sum(wr)
    se = math.sqrt(1 / sum(wr))
    z = pooled / se if se > 0 else math.nan
    return {
        "k": k,
        "pooled_effect": pooled,
        "se": se,
        "ci_low": pooled - 1.96 * se,
        "ci_high": pooled + 1.96 * se,
        "p_value": normal_p_from_z(z) if not math.isnan(z) else math.nan,
        "tau2": tau2,
        "I2": max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0,
    }


def linear_slope(x, y):
    pairs = clean_pairs(x, y)
    if len(pairs) < 3:
        return math.nan
    xs, ys = zip(*pairs)
    mx, my = mean(xs), mean(ys)
    den = sum((a - mx) ** 2 for a in xs)
    if den == 0:
        return math.nan
    return sum((a - mx) * (b - my) for a, b in pairs) / den


def residualize(y, x):
    pairs = clean_pairs(x, y)
    if len(pairs) < 3:
        return [math.nan for _ in y]
    xs, ys = zip(*pairs)
    b = linear_slope(xs, ys)
    a = mean(ys) - b * mean(xs)
    return [yi - (a + b * xi) if not math.isnan(yi) and not math.isnan(xi) else math.nan for xi, yi in zip(x, y)]


def mediation_indirect(x, m, y):
    # a: M ~ X. b: Y residualized on X vs M residualized on X.
    a = linear_slope(x, m)
    y_res = residualize(y, x)
    m_res = residualize(m, x)
    b = linear_slope(m_res, y_res)
    return a, b, a * b if not math.isnan(a) and not math.isnan(b) else math.nan


def bootstrap_ci_indirect(rows, x_col, m_col, y_col, n_boot=1000):
    data = [r for r in rows if not any(math.isnan(fnum(r.get(c))) for c in [x_col, m_col, y_col])]
    n = len(data)
    if n < 30:
        return math.nan, math.nan, math.nan, n
    x = [fnum(r[x_col]) for r in data]
    m = [fnum(r[m_col]) for r in data]
    y = [fnum(r[y_col]) for r in data]
    _, _, point = mediation_indirect(x, m, y)
    boots = []
    for _ in range(n_boot):
        sample = [data[random.randrange(n)] for _ in range(n)]
        xb = [fnum(r[x_col]) for r in sample]
        mb = [fnum(r[m_col]) for r in sample]
        yb = [fnum(r[y_col]) for r in sample]
        _, _, indirect = mediation_indirect(xb, mb, yb)
        if not math.isnan(indirect):
            boots.append(indirect)
    if len(boots) < 50:
        return point, math.nan, math.nan, n
    boots.sort()
    lo = boots[int(0.025 * (len(boots) - 1))]
    hi = boots[int(0.975 * (len(boots) - 1))]
    return point, lo, hi, n


def fmt(x, digits=4):
    if isinstance(x, str):
        return x
    if x is None or math.isnan(x):
        return "NA"
    if abs(x) >= 1000 or (abs(x) < 0.001 and x != 0):
        return f"{x:.3e}"
    return f"{x:.{digits}f}"


def write_docx(path, title, paragraphs):
    content = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>',
        f'<w:p><w:r><w:rPr><w:b/></w:rPr><w:t>{escape(title)}</w:t></w:r></w:p>',
    ]
    for para in paragraphs:
        content.append(f'<w:p><w:r><w:t>{escape(para)}</w:t></w:r></w:p>')
    content.append('<w:sectPr/></w:body></w:document>')
    types = """<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>"""
    rels = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("word/document.xml", "\n".join(content))


with open(LOG_PATH, "w", encoding="utf-8") as log:
    log.write(f"Started {datetime.now().isoformat()}\n")
    log.write(f"Random seed: {RANDOM_SEED}\n")

bulk_summary = read_csv(os.path.join(SRC_BULK, "bulk_meta_interaction_summary.csv"))
bulk_scores = read_csv(os.path.join(SRC_BULK, "bulk_signature_scores_merged.csv"))
cnv_candidates = read_csv(os.path.join(SRC_CNV, "CNV_expression_dosage_genes_positive_rho_padj005.csv"))
cnv_all = read_csv(os.path.join(SRC_CNV, "CNV_expression_dosage_correlation_by_infercnv_subcluster.csv"))

# P2-1 bulk endpoint repositioning: keep negative Cox; add covariation/meta.
bulk_negative_rows = []
for row in bulk_summary:
    beta = fnum(row.get("beta"))
    se = fnum(row.get("se"))
    p = fnum(row.get("p_value"))
    bulk_negative_rows.append({
        "endpoint": row.get("endpoint", ""),
        "model": row.get("model", ""),
        "n_cohorts": row.get("n_cohorts", ""),
        "beta": beta,
        "se": se,
        "HR": row.get("HR", ""),
        "CI_low": row.get("CI_low", ""),
        "CI_high": row.get("CI_high", ""),
        "p_value": p,
        "interpretation": "negative_or_not_significant; OS/PFS not used as primary endpoint",
    })
bh_adjust(bulk_negative_rows)
write_csv(os.path.join(DIRS["bulk"], "bulk_survival_interaction_negative_results_retained.csv"), bulk_negative_rows,
          ["endpoint", "model", "n_cohorts", "beta", "se", "HR", "CI_low", "CI_high", "p_value", "p_adj_BH", "interpretation"])

bulk_corr_pairs = [
    ("SPP1_TAM_score", "ITGB1_CD44_tumor_score"),
    ("SPP1_TAM_score", "KRAS_Hypoxia_score"),
    ("macrophage_fraction", "ITGB1_CD44_tumor_score"),
    ("macrophage_fraction", "KRAS_Hypoxia_score"),
    ("ITGB1_CD44_tumor_score", "KRAS_Hypoxia_score"),
]
cohort_rows = defaultdict(list)
for row in bulk_scores:
    cohort_rows[row.get("cohort", "unknown")].append(row)

bulk_cov_rows = []
for cohort, rows in sorted(cohort_rows.items()):
    for xcol, ycol in bulk_corr_pairs:
        x = [fnum(r.get(xcol)) for r in rows]
        y = [fnum(r.get(ycol)) for r in rows]
        r, p, n = spearman(x, y)
        if n > 3 and not math.isnan(r) and abs(r) < 1:
            fisher_z = 0.5 * math.log((1 + r) / (1 - r))
            variance = 1 / (n - 3)
        else:
            fisher_z = math.nan
            variance = math.nan
        bulk_cov_rows.append({
            "cohort": cohort,
            "x": xcol,
            "y": ycol,
            "n_samples": n,
            "spearman_r": r,
            "p_value": p,
            "fisher_z": fisher_z,
            "variance": variance,
        })
bh_adjust(bulk_cov_rows)
write_csv(os.path.join(DIRS["bulk"], "bulk_signature_covariation_by_cohort.csv"), bulk_cov_rows,
          ["cohort", "x", "y", "n_samples", "spearman_r", "p_value", "p_adj_BH", "fisher_z", "variance"])

bulk_meta_rows = []
for xcol, ycol in bulk_corr_pairs:
    rows = [r for r in bulk_cov_rows if r["x"] == xcol and r["y"] == ycol]
    res = random_effects_meta([fnum(r["fisher_z"]) for r in rows], [fnum(r["variance"]) for r in rows])
    pooled_z = res.get("pooled_effect", math.nan)
    bulk_meta_rows.append({
        "x": xcol,
        "y": ycol,
        "k_cohorts": res.get("k", 0),
        "pooled_spearman_r": math.tanh(pooled_z) if not math.isnan(pooled_z) else math.nan,
        "ci_low_r": math.tanh(res.get("ci_low", math.nan)) if res.get("k", 0) else math.nan,
        "ci_high_r": math.tanh(res.get("ci_high", math.nan)) if res.get("k", 0) else math.nan,
        "p_value": res.get("p_value", math.nan),
        "tau2": res.get("tau2", math.nan),
        "I2": res.get("I2", math.nan),
        "method": "random-effects meta-analysis of cohort-level Spearman Fisher z",
    })
bh_adjust(bulk_meta_rows)
write_csv(os.path.join(DIRS["bulk"], "bulk_signature_covariation_random_effects_meta.csv"), bulk_meta_rows,
          ["x", "y", "k_cohorts", "pooled_spearman_r", "ci_low_r", "ci_high_r", "p_value", "p_adj_BH", "tau2", "I2", "method"])

loo_rows = []
for xcol, ycol in bulk_corr_pairs:
    rows = [r for r in bulk_cov_rows if r["x"] == xcol and r["y"] == ycol]
    full = next((r for r in bulk_meta_rows if r["x"] == xcol and r["y"] == ycol), {})
    for left in sorted(set(r["cohort"] for r in rows)):
        kept = [r for r in rows if r["cohort"] != left]
        res = random_effects_meta([fnum(r["fisher_z"]) for r in kept], [fnum(r["variance"]) for r in kept])
        pooled = math.tanh(res.get("pooled_effect", math.nan)) if res.get("k", 0) else math.nan
        full_r = fnum(full.get("pooled_spearman_r"))
        loo_rows.append({
            "x": xcol,
            "y": ycol,
            "full_pooled_spearman_r": full_r,
            "left_out_cohort": left,
            "loo_k_cohorts": res.get("k", 0),
            "loo_pooled_spearman_r": pooled,
            "direction_preserved": "" if math.isnan(full_r) or math.isnan(pooled) else (full_r == 0 or full_r * pooled > 0),
        })
write_csv(os.path.join(DIRS["bulk"], "bulk_signature_covariation_leave_one_cohort_out.csv"), loo_rows,
          ["x", "y", "full_pooled_spearman_r", "left_out_cohort", "loo_k_cohorts", "loo_pooled_spearman_r", "direction_preserved"])

# P2-2 CNV-expression dosage independent validation status.
renamed = []
for row in cnv_candidates:
    renamed.append({
        "gene": row.get("gene", ""),
        "rho": row.get("rho", ""),
        "p_value": row.get("p_value", ""),
        "padj": row.get("padj", ""),
        "n_groups": row.get("n_groups", ""),
        "cnv_sd": row.get("cnv_sd", ""),
        "expr_sd": row.get("expr_sd", ""),
        "evidence_label": "RNA-derived CNV-expression coupling candidate",
        "independent_DNA_RNA_validation_status": "not_validated_in_local_outputs",
        "interpretation": "do not call DNA dosage driver without independent DNA CNV plus RNA validation",
    })
write_csv(os.path.join(DIRS["evidence"], "RNA_derived_CNV_expression_coupling_candidates.csv"), renamed,
          ["gene", "rho", "p_value", "padj", "n_groups", "cnv_sd", "expr_sd", "evidence_label",
           "independent_DNA_RNA_validation_status", "interpretation"])

cnv_summary = [{
    "source_table": "CNV_expression_dosage_genes_positive_rho_padj005.csv",
    "n_positive_candidates_padj005": len(renamed),
    "top_gene": renamed[0]["gene"] if renamed else "",
    "top_rho": renamed[0]["rho"] if renamed else "",
    "top_padj": renamed[0]["padj"] if renamed else "",
    "required_for_upgrade": "TCGA/GISTIC SNP-array or paired WES/WGS CNV plus matched RNA validation",
    "current_label": "RNA-derived CNV-expression coupling candidates",
}]
write_csv(os.path.join(DIRS["evidence"], "cnv_expression_independent_validation_status.csv"), cnv_summary,
          ["source_table", "n_positive_candidates_padj005", "top_gene", "top_rho", "top_padj",
           "required_for_upgrade", "current_label"])

# P2-3 exploratory mediation in bulk only.
med_rows = []
for cohort, rows in sorted(cohort_rows.items()):
    point, lo, hi, n = bootstrap_ci_indirect(rows, "SPP1_TAM_score", "ITGB1_CD44_tumor_score", "KRAS_Hypoxia_score", n_boot=1000)
    med_rows.append({
        "cohort": cohort,
        "n_samples": n,
        "x": "SPP1_TAM_score",
        "mediator": "ITGB1_CD44_tumor_score",
        "outcome": "KRAS_Hypoxia_score",
        "indirect_effect": point,
        "bootstrap_ci_low": lo,
        "bootstrap_ci_high": hi,
        "status": "exploratory_bulk_statistical_mediation" if n >= 30 else "not_run_insufficient_samples",
        "interpretation": "statistical mediation only; not causal and not spatial niche proof",
    })
write_csv(os.path.join(DIRS["evidence"], "bulk_exploratory_mediation_bootstrap.csv"), med_rows,
          ["cohort", "n_samples", "x", "mediator", "outcome", "indirect_effect", "bootstrap_ci_low", "bootstrap_ci_high",
           "status", "interpretation"])

# Integrated evidence scoring.
p1_status_path = os.path.join(DIRS["tables"], "p1_execution_status.csv")
p1_status = read_csv(p1_status_path)
spatial_meta = read_csv(os.path.join(OUT_ROOT, "09_spatial_deconvolution_neighborhood", "spatial_random_effects_meta_analysis.csv"))
external_assoc = read_csv(os.path.join(OUT_ROOT, "06_external_scrna_projection", "external_scrna_patient_level_associations.csv"))

bulk_sp = next((r for r in bulk_meta_rows if r["x"] == "SPP1_TAM_score" and r["y"] == "KRAS_Hypoxia_score"), {})
bulk_med_sig = [r for r in med_rows if fnum(r.get("bootstrap_ci_low")) > 0 or fnum(r.get("bootstrap_ci_high")) < 0]
external_0204 = next((r for r in external_assoc if r.get("stratum") == "all_samples" and r.get("y") == "Subclone02_04_common_mean_target_tumor" and r.get("x") == "source_myeloid_fraction"), {})
spatial_all = next((r for r in spatial_meta if r.get("dataset") == "ALL"), {})

evidence_rows = [
    {
        "evidence_layer": "Gate A clone cross-patient validity",
        "classification": "因数据结构无法判断",
        "key_stat": "integrated_oc patient_id missing",
        "limitation": "cannot call old CNV_Subclone labels cross-patient states",
    },
    {
        "evidence_layer": "External SPP1-myeloid vs 02/04-like association",
        "classification": "未支持",
        "key_stat": f"Spearman r={fmt(fnum(external_0204.get('spearman_r')))}, BH q={fmt(fnum(external_0204.get('p_adj_BH')))}",
        "limitation": "external signature projection; site/treatment may be confounded",
    },
    {
        "evidence_layer": "Spatial co-expression meta",
        "classification": "部分支持" if fnum(spatial_all.get("p_adj_BH")) < 0.05 else "未支持",
        "key_stat": f"pooled r={fmt(fnum(spatial_all.get('pooled_spearman_r')))}, BH q={fmt(fnum(spatial_all.get('p_adj_BH')))}",
        "limitation": "high heterogeneity; not full spatial deconvolution",
    },
    {
        "evidence_layer": "Bulk OS/PFS interaction",
        "classification": "未支持",
        "key_stat": "Cox interaction models not significant",
        "limitation": "OS/PFS not primary endpoint",
    },
    {
        "evidence_layer": "Bulk SPP1/KRAS-hypoxia covariation",
        "classification": "部分支持" if fnum(bulk_sp.get("p_adj_BH")) < 0.05 else "未支持",
        "key_stat": f"pooled r={fmt(fnum(bulk_sp.get('pooled_spearman_r')))}, BH q={fmt(fnum(bulk_sp.get('p_adj_BH')))}",
        "limitation": "bulk signatures are not local cell-cell interaction evidence",
    },
    {
        "evidence_layer": "CNV-expression dosage",
        "classification": "部分支持",
        "key_stat": f"{len(renamed)} RNA-derived candidates; top={cnv_summary[0]['top_gene']} rho={cnv_summary[0]['top_rho']}",
        "limitation": "RNA-derived inferCNV circularity; no independent DNA CNV validation found locally",
    },
    {
        "evidence_layer": "Exploratory mediation",
        "classification": "部分支持" if bulk_med_sig else "未支持",
        "key_stat": f"{len(bulk_med_sig)} cohorts with bootstrap CI excluding 0",
        "limitation": "statistical mediation only; not causal",
    },
    {
        "evidence_layer": "NicheNet mechanism",
        "classification": "未支持",
        "key_stat": "full NicheNet output absent",
        "limitation": "requires patient-level stable 02/04 DEG target set",
    },
]
write_csv(os.path.join(DIRS["evidence"], "integrated_evidence_scoring.csv"), evidence_rows,
          ["evidence_layer", "classification", "key_stat", "limitation"])

bulk_report = f"""# P2 bulk, CNV-expression and mediation report

Generated: {datetime.now().isoformat(timespec="seconds")}

## P2-1 Bulk endpoint repositioning

The five-cohort Cox interaction results were retained as negative or non-significant and are not used as the primary endpoint.

Output: `11_bulk_negative_validation/bulk_survival_interaction_negative_results_retained.csv`

Bulk signature covariation was re-positioned as a co-variation check. Main meta result for `SPP1_TAM_score` vs `KRAS_Hypoxia_score`:

- pooled Spearman r = {fmt(fnum(bulk_sp.get('pooled_spearman_r')))}
- 95% CI = [{fmt(fnum(bulk_sp.get('ci_low_r')))}, {fmt(fnum(bulk_sp.get('ci_high_r')))}]
- BH q = {fmt(fnum(bulk_sp.get('p_adj_BH')))}
- I2 = {fmt(fnum(bulk_sp.get('I2')))}

Output: `11_bulk_negative_validation/bulk_signature_covariation_random_effects_meta.csv`

## P2-2 CNV-expression dosage independent validation status

Local outputs contain {len(renamed)} positive RNA-derived CNV-expression coupling candidates at the existing threshold. These were relabelled as:

`RNA-derived CNV-expression coupling candidates`

They were not upgraded to DNA dosage drivers because no independent DNA CNV plus RNA validation table was found locally.

Output: `12_integrated_evidence_scoring/RNA_derived_CNV_expression_coupling_candidates.csv`

## P2-3 Statistical mediation exploration

Bulk-level exploratory mediation was run per cohort:

`SPP1_TAM_score -> ITGB1_CD44_tumor_score -> KRAS_Hypoxia_score`

This is explicitly statistical mediation, not causal proof and not local spatial niche evidence. Cohorts with bootstrap CI excluding 0: {len(bulk_med_sig)}.

Output: `12_integrated_evidence_scoring/bulk_exploratory_mediation_bootstrap.csv`

## Integrated conclusion

The current evidence supports cautious, mixed candidate biology rather than a validated mechanism. The strongest limitation remains unresolved P0/Gate A: integrated_oc lacks `patient_id`, so current CNV_Subclone_02/04 cannot yet be promoted to cross-patient CNV programs.
"""

with open(os.path.join(DIRS["reports"], "05_integrated_evidence_and_limitations.md"), "w", encoding="utf-8") as handle:
    handle.write(bulk_report)

html = bulk_report.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br>\n")
with open(os.path.join(DIRS["reports"], "05_integrated_evidence_and_limitations.html"), "w", encoding="utf-8") as handle:
    handle.write("<!doctype html><meta charset='utf-8'><title>P2 integrated evidence</title>"
                 "<body style='font-family:Arial,Microsoft YaHei,sans-serif;line-height:1.55;margin:32px'>"
                 f"{html}</body>")

docx_paragraphs = [
    "P2 bulk endpoint results were retained as negative/non-significant; OS/PFS is not used as the primary endpoint.",
    f"Bulk SPP1_TAM_score vs KRAS_Hypoxia_score meta result: pooled Spearman r = {fmt(fnum(bulk_sp.get('pooled_spearman_r')))}, BH q = {fmt(fnum(bulk_sp.get('p_adj_BH')))}.",
    f"{len(renamed)} RNA-derived CNV-expression coupling candidates were relabelled and not called DNA dosage drivers.",
    f"Exploratory bulk mediation found {len(bulk_med_sig)} cohorts with bootstrap CI excluding 0; this is statistical mediation only, not causal proof.",
    "Integrated conclusion: evidence is mixed/partial; Gate A remains unresolved because patient_id is missing in integrated_oc metadata.",
]
write_docx(os.path.join(DIRS["reports"], "05_integrated_evidence_and_limitations.docx"),
           "P2 Integrated Evidence and Limitations", docx_paragraphs)

summary_path = os.path.join(DIRS["reports"], "FINAL_CODEX_EXECUTION_SUMMARY.md")
summary_text = ""
if os.path.exists(summary_path):
    with open(summary_path, "r", encoding="utf-8") as handle:
        summary_text = handle.read()
marker = "\n## P2 execution update\n"
if marker in summary_text:
    summary_text = summary_text.split(marker)[0].rstrip() + "\n"
with open(summary_path, "w", encoding="utf-8") as handle:
    handle.write(summary_text.rstrip())
    handle.write("\n\n## P2 execution update\n\n")
    handle.write(bulk_report)

with open(LOG_PATH, "a", encoding="utf-8") as log:
    log.write(f"bulk summary rows: {len(bulk_summary)}\n")
    log.write(f"bulk score rows: {len(bulk_scores)}\n")
    log.write(f"cnv candidate rows: {len(renamed)}\n")
    log.write(f"finished: {datetime.now().isoformat()}\n")
