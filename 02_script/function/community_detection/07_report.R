# ============================================================
# 07_report.R
# Génération d'un rapport HTML interactif pour les top candidats
#
# Produit un fichier HTML auto-contenu (images en base64, CSS/JS inline)
# qui permet de comparer visuellement et interactivement les finalistes
# du grid search avant de choisir le clustering final.
#
# Usage :
#   source("07_report.R")
#   generate_report(
#     all_results = results,          # data.frame issu de run_grid_search()
#     base_dir    = "path/to/outputs",
#     out_html    = "clustering_report.html",
#     n_top       = 5
#   )
# ============================================================

# ── Dépendances ───────────────────────────────────────────────
if (!requireNamespace("base64enc", quietly = TRUE)) install.packages("base64enc")
library(base64enc)

# ── Encode une image PNG en data URI (pour l'inliner dans le HTML) ────────────
img_to_base64 <- function(path) {
  if (!file.exists(path)) return(NULL)
  paste0("data:image/png;base64,",
         base64enc::base64encode(path))
}

# ── Cherche un fichier PNG dans un sous-dossier (flexible) ───────────────────
find_png <- function(base, ...) {
  p <- file.path(base, ...)
  if (file.exists(p)) return(p)
  NULL
}

# ── Collecte tous les assets d'un combo ──────────────────────────────────────
collect_combo_assets <- function(base_dir, k, dims_config, resolution) {
  combo_dir <- file.path(base_dir, sprintf("k%d_%s", k, dims_config))
  res_dir   <- file.path(combo_dir, sprintf("res_%.2f", resolution))
  bio_dir   <- file.path(res_dir, "Bio")

  list(
    # Métriques (lues depuis le data.frame, passées séparément)
    # Plots mathématiques
    umap_density   = find_png(res_dir, "UMAP", "UMAP_density.png"),
    umap_sil       = find_png(res_dir, "UMAP", "UMAP_silhouette.png"),
    umap_facet     = find_png(res_dir, "UMAP", "UMAP_density_facet.png"),
    sil_violin     = find_png(res_dir, "Silhouette", "silhouette_violin.png"),
    sil_mean       = find_png(res_dir, "Silhouette", "silhouette_mean_per_cluster.png"),
    cluster_sizes  = find_png(res_dir, "Cluster_sizes", "cluster_cell_counts.png"),
    batch          = find_png(res_dir, "Batch", "batch_cluster_proportion.png"),
    qc_violin      = find_png(res_dir, "QC", "QC_violin.png"),
    qc_mito        = find_png(res_dir, "QC", "QC_mito_ridges.png"),
    qc_cc          = find_png(res_dir, "QC", "QC_cell_cycle.png"),
    clustree       = find_png(combo_dir, "Clustree", sprintf("clustree_k%d.png", k)),
    # Plots biologiques (Phase 2 — peuvent être absents)
    heatmap        = find_png(bio_dir, "Heatmap", "top_markers_heatmap.png"),
    # Module scores : liste de tous les PNGs trouvés
    module_scores  = {
      md <- file.path(bio_dir, "ModuleScores")
      if (dir.exists(md)) list.files(md, pattern = "\\.png$", full.names = TRUE)
      else character(0)
    }
  )
}

# ── Génère le bloc HTML d'un candidat ────────────────────────────────────────
render_candidate_card <- function(row, assets, rank) {
  label <- sprintf("k=%d | %s | res=%.2f", row$k, row$dims_config, row$resolution)

  # Score badges
  score_html <- function(label, value, color) {
    if (is.na(value)) return("")
    sprintf('<div class="badge" style="--badge-color:%s"><span class="badge-label">%s</span><span class="badge-value">%.3f</span></div>',
            color, label, value)
  }

  badges <- paste0(
    score_html("Silhouette",  row$mean_sil,      "#4a9eff"),
    score_html("ARI",         row$mean_ari,       "#7ec8a0"),
    score_html("Math score",  row$math_score,     "#f0a500"),
    score_html("Bio score",   row$bio_score,      "#e07070"),
    score_html("Final score", row$final_score,    "#c084fc")
  )

  meta_html <- sprintf(
    '<div class="meta-grid">
      <div class="meta-item"><span class="meta-key">k</span><span class="meta-val">%d</span></div>
      <div class="meta-item"><span class="meta-key">dims</span><span class="meta-val">%s</span></div>
      <div class="meta-item"><span class="meta-key">res</span><span class="meta-val">%.2f</span></div>
      <div class="meta-item"><span class="meta-key">clusters</span><span class="meta-val">%d</span></div>
      <div class="meta-item"><span class="meta-key">prop_neg</span><span class="meta-val">%.1f%%</span></div>
      <div class="meta-item"><span class="meta-key">sd_ARI</span><span class="meta-val">%.3f</span></div>
    </div>',
    row$k, row$dims_config, row$resolution, row$n_clusters,
    100 * row$prop_negative,
    ifelse(is.na(row$sd_ari), 0, row$sd_ari)
  )

  # Encode images
  img_tag <- function(path, caption, cls = "plot-img") {
    if (is.null(path) || !file.exists(path)) return("")
    b64 <- img_to_base64(path)
    if (is.null(b64)) return("")
    sprintf('<figure class="plot-figure %s">
               <img src="%s" alt="%s" loading="lazy">
               <figcaption>%s</figcaption>
             </figure>', cls, b64, caption, caption)
  }

  # Module score thumbnails
  module_imgs <- if (length(assets$module_scores) > 0) {
    paste(sapply(assets$module_scores, function(p) {
      sig_name <- gsub("ModuleScore_UMAP_ModuleScore_|ModuleScore_UMAP_|\\.png$", "", basename(p))
      img_tag(p, sig_name, "plot-img plot-thumb")
    }), collapse = "\n")
  } else ""

  sprintf('
  <div class="candidate-card" id="cand-%d" data-rank="%d" data-final-score="%.4f">
    <div class="card-header">
      <div class="rank-badge">%s%d</div>
      <h2 class="card-title">%s</h2>
    </div>

    <div class="scores-row">%s</div>
    %s

    <div class="tabs" data-card="%d">
      <button class="tab active" data-tab="umap">UMAP</button>
      <button class="tab" data-tab="silhouette">Silhouette</button>
      <button class="tab" data-tab="qc">QC</button>
      <button class="tab" data-tab="clustree">Clustree</button>
      <button class="tab" data-tab="bio">Biologie</button>
    </div>

    <div class="tab-content active" data-card="%d" data-panel="umap">
      <div class="plot-grid">
        %s %s %s
      </div>
    </div>

    <div class="tab-content" data-card="%d" data-panel="silhouette">
      <div class="plot-grid">
        %s %s
      </div>
    </div>

    <div class="tab-content" data-card="%d" data-panel="qc">
      <div class="plot-grid">
        %s %s %s %s
      </div>
    </div>

    <div class="tab-content" data-card="%d" data-panel="clustree">
      <div class="plot-grid">
        %s
      </div>
    </div>

    <div class="tab-content" data-card="%d" data-panel="bio">
      <div class="plot-grid">
        %s
        %s
        <div class="module-scores-grid">%s</div>
      </div>
    </div>

  </div>',
    rank, rank, ifelse(is.na(row$final_score), 0, row$final_score),
    ifelse(rank == 1, "★ ", ""), rank,
    label,
    badges,
    meta_html,
    rank, rank,
    # UMAP tab
    img_tag(assets$umap_density, "UMAP — clusters + density"),
    img_tag(assets$umap_sil,     "UMAP — silhouette score"),
    img_tag(assets$umap_facet,   "UMAP — par condition"),
    # Silhouette tab
    rank,
    img_tag(assets$sil_violin, "Silhouette par cluster"),
    img_tag(assets$sil_mean,   "Silhouette moyenne par cluster"),
    # QC tab
    rank,
    img_tag(assets$qc_violin,  "QC metrics (nFeature, nCount, %mt)"),
    img_tag(assets$qc_mito,    "% Mitochondrial — ridges"),
    img_tag(assets$qc_cc,      "Cycle cellulaire par cluster"),
    img_tag(assets$cluster_sizes, "Taille des clusters"),
    # Clustree tab
    rank,
    img_tag(assets$clustree, "Clustree — stabilité des clusters"),
    # Bio tab
    rank,
    img_tag(assets$heatmap,  "Heatmap top markers"),
    "",  # placeholder volcano (trop de fichiers à lister)
    module_imgs
  )
}

# ── Template HTML complet ─────────────────────────────────────────────────────
html_template <- function(candidates_html, summary_table_html,
                          score_radar_json, n_top, timestamp) {
  sprintf('<!DOCTYPE html>
<html lang="fr">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Clustering Report — Top %d candidats</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:ital,wght@0,300;0,500;0,700;1,300&display=swap" rel="stylesheet">
<style>
:root {
  --bg:        #0d0f14;
  --surface:   #13161e;
  --surface2:  #1a1e2a;
  --border:    #252a38;
  --accent:    #4a9eff;
  --accent2:   #c084fc;
  --gold:      #f0a500;
  --green:     #7ec8a0;
  --red:       #e07070;
  --text:      #e8eaf0;
  --muted:     #6b7280;
  --mono:      "Space Mono", monospace;
  --sans:      "DM Sans", sans-serif;
  --radius:    12px;
  --transition: 200ms cubic-bezier(0.4, 0, 0.2, 1);
}

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  background: var(--bg);
  color: var(--text);
  font-family: var(--sans);
  font-size: 14px;
  line-height: 1.6;
  min-height: 100vh;
}

/* ── Header ── */
.site-header {
  background: linear-gradient(135deg, #0d0f14 0%%, #13161e 60%%, #1a1432 100%%);
  border-bottom: 1px solid var(--border);
  padding: 2.5rem 3rem 2rem;
  position: sticky;
  top: 0;
  z-index: 100;
  backdrop-filter: blur(12px);
}

.header-top {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: 2rem;
  flex-wrap: wrap;
}

.header-eyebrow {
  font-family: var(--mono);
  font-size: 10px;
  letter-spacing: 0.15em;
  color: var(--accent);
  text-transform: uppercase;
  margin-bottom: 0.4rem;
}

h1 {
  font-family: var(--mono);
  font-size: clamp(1.2rem, 3vw, 1.8rem);
  font-weight: 700;
  color: var(--text);
  letter-spacing: -0.02em;
}

.header-meta {
  font-size: 12px;
  color: var(--muted);
  font-family: var(--mono);
  margin-top: 0.4rem;
}

/* ── Navigation compacte ── */
.nav-pills {
  display: flex;
  gap: 0.5rem;
  flex-wrap: wrap;
  margin-top: 1.5rem;
}

.nav-pill {
  background: var(--surface2);
  border: 1px solid var(--border);
  color: var(--muted);
  font-family: var(--mono);
  font-size: 11px;
  padding: 0.3rem 0.8rem;
  border-radius: 100px;
  cursor: pointer;
  transition: var(--transition);
  text-decoration: none;
}

.nav-pill:hover, .nav-pill.active {
  background: var(--accent);
  border-color: var(--accent);
  color: #fff;
}

/* ── Layout ── */
.main {
  max-width: 1600px;
  margin: 0 auto;
  padding: 2rem 2rem 4rem;
}

/* ── Summary Table ── */
.section-title {
  font-family: var(--mono);
  font-size: 11px;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--muted);
  margin-bottom: 1rem;
  padding-bottom: 0.5rem;
  border-bottom: 1px solid var(--border);
}

.summary-section {
  margin-bottom: 3rem;
}

.summary-table {
  width: 100%%;
  border-collapse: collapse;
  font-family: var(--mono);
  font-size: 12px;
}

.summary-table th {
  background: var(--surface2);
  color: var(--muted);
  font-size: 10px;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  padding: 0.6rem 0.9rem;
  text-align: left;
  border-bottom: 1px solid var(--border);
}

.summary-table td {
  padding: 0.55rem 0.9rem;
  border-bottom: 1px solid var(--border);
  color: var(--text);
}

.summary-table tr:hover td { background: var(--surface2); }
.summary-table tr.top-row td { border-left: 3px solid var(--gold); }

.score-bar {
  display: inline-flex;
  align-items: center;
  gap: 0.5rem;
  width: 100%%;
}

.score-bar-inner {
  height: 4px;
  border-radius: 2px;
  background: var(--accent);
  transition: width 0.6s ease;
}

/* ── Candidate cards ── */
.candidates-grid {
  display: flex;
  flex-direction: column;
  gap: 2.5rem;
}

.candidate-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  overflow: hidden;
  transition: var(--transition);
  scroll-margin-top: 120px;
}

.candidate-card:hover {
  border-color: #2d3452;
}

.card-header {
  display: flex;
  align-items: center;
  gap: 1.2rem;
  padding: 1.4rem 1.6rem 1rem;
  background: linear-gradient(90deg, var(--surface2) 0%%, transparent 100%%);
  border-bottom: 1px solid var(--border);
}

.rank-badge {
  font-family: var(--mono);
  font-size: 13px;
  font-weight: 700;
  background: var(--surface);
  border: 1px solid var(--border);
  color: var(--accent);
  width: 36px;
  height: 36px;
  border-radius: 50%%;
  display: flex;
  align-items: center;
  justify-content: center;
  flex-shrink: 0;
}

.candidate-card[data-rank="1"] .rank-badge {
  background: var(--gold);
  border-color: var(--gold);
  color: #000;
}

.card-title {
  font-family: var(--mono);
  font-size: 15px;
  font-weight: 700;
  color: var(--text);
  letter-spacing: -0.01em;
}

/* ── Score badges ── */
.scores-row {
  display: flex;
  flex-wrap: wrap;
  gap: 0.5rem;
  padding: 1rem 1.6rem;
  border-bottom: 1px solid var(--border);
}

.badge {
  display: flex;
  flex-direction: column;
  align-items: center;
  background: color-mix(in srgb, var(--badge-color) 10%%, var(--surface2));
  border: 1px solid color-mix(in srgb, var(--badge-color) 30%%, transparent);
  border-radius: 8px;
  padding: 0.4rem 0.8rem;
  min-width: 80px;
}

.badge-label {
  font-family: var(--mono);
  font-size: 9px;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  color: var(--badge-color);
  opacity: 0.8;
}

.badge-value {
  font-family: var(--mono);
  font-size: 16px;
  font-weight: 700;
  color: var(--badge-color);
  margin-top: 1px;
}

/* ── Meta grid ── */
.meta-grid {
  display: flex;
  flex-wrap: wrap;
  gap: 0;
  padding: 0.8rem 1.6rem;
  border-bottom: 1px solid var(--border);
  background: var(--surface2);
}

.meta-item {
  display: flex;
  flex-direction: column;
  padding: 0.3rem 1.2rem;
  border-right: 1px solid var(--border);
}

.meta-item:last-child { border-right: none; }

.meta-key {
  font-family: var(--mono);
  font-size: 9px;
  text-transform: uppercase;
  letter-spacing: 0.1em;
  color: var(--muted);
}

.meta-val {
  font-family: var(--mono);
  font-size: 15px;
  font-weight: 700;
  color: var(--text);
}

/* ── Tabs ── */
.tabs {
  display: flex;
  gap: 0;
  border-bottom: 1px solid var(--border);
  padding: 0 1.2rem;
  background: var(--surface2);
  overflow-x: auto;
}

.tab {
  background: none;
  border: none;
  border-bottom: 2px solid transparent;
  color: var(--muted);
  font-family: var(--mono);
  font-size: 11px;
  letter-spacing: 0.05em;
  text-transform: uppercase;
  padding: 0.75rem 1rem;
  cursor: pointer;
  transition: var(--transition);
  white-space: nowrap;
}

.tab:hover { color: var(--text); }

.tab.active {
  color: var(--accent);
  border-bottom-color: var(--accent);
}

/* ── Tab content ── */
.tab-content {
  display: none;
  padding: 1.5rem;
}

.tab-content.active { display: block; }

/* ── Plot grid ── */
.plot-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
  gap: 1rem;
}

.plot-figure {
  background: var(--surface2);
  border: 1px solid var(--border);
  border-radius: 8px;
  overflow: hidden;
}

.plot-figure img {
  width: 100%%;
  height: auto;
  display: block;
  cursor: zoom-in;
  transition: opacity var(--transition);
}

.plot-figure img:hover { opacity: 0.9; }

figcaption {
  font-family: var(--mono);
  font-size: 10px;
  color: var(--muted);
  padding: 0.4rem 0.7rem;
  letter-spacing: 0.05em;
  text-transform: uppercase;
  border-top: 1px solid var(--border);
}

.plot-thumb img { height: 200px; object-fit: cover; }

.module-scores-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 0.75rem;
  grid-column: 1 / -1;
}

/* ── Lightbox ── */
.lightbox {
  display: none;
  position: fixed;
  inset: 0;
  background: rgba(0,0,0,0.92);
  z-index: 1000;
  align-items: center;
  justify-content: center;
  cursor: zoom-out;
}

.lightbox.open { display: flex; }

.lightbox img {
  max-width: 92vw;
  max-height: 92vh;
  object-fit: contain;
  border-radius: 4px;
  box-shadow: 0 25px 80px rgba(0,0,0,0.8);
}

/* ── Floating compare bar ── */
.compare-bar {
  position: fixed;
  bottom: 2rem;
  left: 50%%;
  transform: translateX(-50%%);
  background: var(--surface2);
  border: 1px solid var(--border);
  border-radius: 100px;
  padding: 0.6rem 1.4rem;
  display: flex;
  align-items: center;
  gap: 1rem;
  font-family: var(--mono);
  font-size: 11px;
  color: var(--muted);
  box-shadow: 0 8px 32px rgba(0,0,0,0.5);
  opacity: 0;
  pointer-events: none;
  transition: opacity 0.3s;
  z-index: 200;
}

.compare-bar.visible {
  opacity: 1;
  pointer-events: auto;
}

/* ── Scrollbar ── */
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: var(--bg); }
::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }

/* ── Responsive ── */
@media (max-width: 768px) {
  .site-header { padding: 1.5rem; }
  .main { padding: 1rem; }
  .plot-grid { grid-template-columns: 1fr; }
  .meta-grid { gap: 0.5rem; }
}
</style>
</head>
<body>

<header class="site-header">
  <div class="header-top">
    <div>
      <div class="header-eyebrow">scRNA-seq / CITE-seq — WNN Leiden Clustering</div>
      <h1>Clustering Report — Top %d candidats</h1>
      <div class="header-meta">Généré le %s &nbsp;|&nbsp; Grid search terminé</div>
    </div>
  </div>
  <nav class="nav-pills" id="nav-pills">
    <a href="#summary" class="nav-pill active">Vue d\'ensemble</a>
    %s
  </nav>
</header>

<main class="main">

  <!-- Summary table -->
  <section class="summary-section" id="summary">
    <div class="section-title">Tableau de bord — Tous les candidats</div>
    %s
  </section>

  <!-- Candidate cards -->
  <section id="candidates">
    <div class="section-title">Détail par candidat</div>
    <div class="candidates-grid">
      %s
    </div>
  </section>

</main>

<!-- Lightbox -->
<div class="lightbox" id="lightbox">
  <img id="lightbox-img" src="" alt="">
</div>

<!-- Floating hint -->
<div class="compare-bar" id="compare-bar">
  <span>💡</span>
  <span id="compare-msg">Cliquez sur une image pour zoomer</span>
</div>

<script>
// ── Tabs ──────────────────────────────────────────────────────────────────
document.querySelectorAll(".tabs").forEach(tabGroup => {
  const card = tabGroup.dataset.card;
  tabGroup.querySelectorAll(".tab").forEach(tab => {
    tab.addEventListener("click", () => {
      tabGroup.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
      tab.classList.add("active");
      const panel = tab.dataset.tab;
      document.querySelectorAll(`.tab-content[data-card="${card}"]`).forEach(tc => {
        tc.classList.toggle("active", tc.dataset.panel === panel);
      });
    });
  });
});

// ── Nav pills scroll ──────────────────────────────────────────────────────
document.querySelectorAll(".nav-pill[data-target]").forEach(pill => {
  pill.addEventListener("click", e => {
    e.preventDefault();
    document.querySelectorAll(".nav-pill").forEach(p => p.classList.remove("active"));
    pill.classList.add("active");
    const target = document.getElementById(pill.dataset.target);
    if (target) target.scrollIntoView({ behavior: "smooth" });
  });
});

// ── Lightbox ──────────────────────────────────────────────────────────────
const lightbox    = document.getElementById("lightbox");
const lightboxImg = document.getElementById("lightbox-img");

document.querySelectorAll(".plot-figure img").forEach(img => {
  img.addEventListener("click", () => {
    lightboxImg.src = img.src;
    lightbox.classList.add("open");
  });
});

lightbox.addEventListener("click", () => lightbox.classList.remove("open"));

document.addEventListener("keydown", e => {
  if (e.key === "Escape") lightbox.classList.remove("open");
});

// ── Score bars animation ──────────────────────────────────────────────────
const obs = new IntersectionObserver(entries => {
  entries.forEach(entry => {
    if (entry.isIntersecting) {
      entry.target.querySelectorAll(".score-bar-inner").forEach(bar => {
        bar.style.width = bar.dataset.width + "%%";
      });
    }
  });
}, { threshold: 0.1 });

document.querySelectorAll(".summary-table tr").forEach(row => obs.observe(row));

// ── Compare bar hint ──────────────────────────────────────────────────────
const compareBar = document.getElementById("compare-bar");
let hintShown = false;

window.addEventListener("scroll", () => {
  if (!hintShown && window.scrollY > 300) {
    compareBar.classList.add("visible");
    hintShown = true;
    setTimeout(() => compareBar.classList.remove("visible"), 4000);
  }
});
</script>
</body>
</html>',
    n_top, n_top, timestamp,
    # nav pills (remplis dynamiquement)
    "%NAV_PILLS%",
    # summary table
    "%SUMMARY_TABLE%",
    # candidates
    "%CANDIDATES%"
  )
}

# ── Génère la table de résumé ─────────────────────────────────────────────────
render_summary_table <- function(results) {
  max_final <- max(results$final_score, na.rm = TRUE)

  rows <- paste(sapply(seq_len(nrow(results)), function(i) {
    row   <- results[i, ]
    is_top <- i <= 5

    bar <- function(val, color) {
      if (is.na(val)) return("—")
      w <- round(100 * val / max(max_final, 0.01))
      sprintf('<div class="score-bar">
                 <div class="score-bar-inner" data-width="%d" style="width:0;background:%s;flex:0 0 %dpx"></div>
                 <span>%.3f</span>
               </div>', w, color, min(w, 80), val)
    }

    sprintf('<tr class="%s">
      <td style="font-family:var(--mono);font-weight:700;color:var(--accent)">%s%d</td>
      <td style="font-family:var(--mono)">%d</td>
      <td style="font-family:var(--mono)">%s</td>
      <td style="font-family:var(--mono)">%.2f</td>
      <td>%d</td>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
    </tr>',
      ifelse(is_top, "top-row", ""),
      ifelse(i == 1, "★ ", ""), i,
      row$k, row$dims_config, row$resolution, row$n_clusters,
      bar(row$mean_sil,    "#4a9eff"),
      bar(row$mean_ari,    "#7ec8a0"),
      bar(row$math_score,  "#f0a500"),
      bar(row$bio_score,   "#e07070"),
      bar(row$final_score, "#c084fc")
    )
  }), collapse = "\n")

  sprintf('<table class="summary-table">
    <thead><tr>
      <th>Rang</th><th>k</th><th>dims</th><th>res</th><th>clusters</th>
      <th>Silhouette</th><th>ARI</th><th>Math score</th>
      <th>Bio score</th><th>Final score</th>
    </tr></thead>
    <tbody>%s</tbody>
  </table>', rows)
}

# ── Point d'entrée principal ──────────────────────────────────────────────────
#
# @param all_results  data.frame issu de run_grid_search(), trié par final_score desc
# @param base_dir     Répertoire racine du grid search (là où sont les plots)
# @param out_html     Chemin du fichier HTML à générer
# @param n_top        Nombre de candidats à inclure dans le rapport (défaut 5)
# ─────────────────────────────────────────────────────────────────────────────
generate_report <- function(all_results, base_dir, out_html, n_top = 5) {
  message(sprintf("[Report] Génération du rapport HTML — top %d candidats...", n_top))

  top <- head(all_results, n_top)

  # Nav pills
  nav_pills_html <- paste(sapply(seq_len(nrow(top)), function(i) {
    row <- top[i, ]
    sprintf('<a href="#cand-%d" class="nav-pill" data-target="cand-%d">%s%d — k%d res%.2f</a>',
            i, i, ifelse(i == 1, "★ ", ""), i, row$k, row$resolution)
  }), collapse = "\n")

  # Summary table
  summary_html <- render_summary_table(all_results)

  # Candidate cards
  cards_html <- paste(sapply(seq_len(nrow(top)), function(i) {
    row    <- top[i, ]
    assets <- collect_combo_assets(base_dir, row$k, row$dims_config, row$resolution)
    message(sprintf("  [Report] Encodage candidat %d/%d (k=%d | res=%.2f)...",
                    i, n_top, row$k, row$resolution))
    render_candidate_card(row, assets, i)
  }), collapse = "\n")

  # Assembly
  timestamp <- format(Sys.time(), "%d/%m/%Y à %H:%M")
  html <- html_template("", "", "", n_top, timestamp)
  html <- gsub("%NAV_PILLS%",    nav_pills_html, html, fixed = TRUE)
  html <- gsub("%SUMMARY_TABLE%", summary_html,  html, fixed = TRUE)
  html <- gsub("%CANDIDATES%",   cards_html,     html, fixed = TRUE)

  writeLines(html, out_html)
  message(sprintf("[Report] ✓ Rapport généré : %s", out_html))
  invisible(out_html)
}

# ── Exemple d'usage ───────────────────────────────────────────────────────────
# generate_report(
#   all_results = results,
#   base_dir    = file.path(DIR_RESULT, "Community_detection", "grid_search"),
#   out_html    = file.path(DIR_RESULT, "Community_detection", "clustering_report.html"),
#   n_top       = 5
# )
