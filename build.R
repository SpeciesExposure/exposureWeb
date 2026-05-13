# =============================================================================
# build.R
#
# Run this script from apps/dashboard/ to prepare data and render the
# dashboard to ../../docs/ (repo root docs/ for GitHub Pages).
#
# Usage (from apps/dashboard/):
#   Rscript build.R
#
# Optional: limit to first N species for fast dev iteration:
#   DEV_N_SPECIES=10 Rscript build.R
# =============================================================================

# Ensure working directory is apps/dashboard/
if (!file.exists("_quarto.yml")) {
  stop("Run this script from apps/dashboard/ (where _quarto.yml lives).")
}

# ---------------------------------------------------------------------------
# Data location: all input and output files live in ../../data/ (repo root data/).
# Set DATA_DIR / CACHE_DIR environment variables to override.
# ---------------------------------------------------------------------------
Sys.setenv(DATA_DIR  = Sys.getenv("DATA_DIR",  "../../data"))
Sys.setenv(CACHE_DIR = Sys.getenv("CACHE_DIR", "data"))
Sys.setenv(DEV_N_SPECIES = Sys.getenv("DEV_N_SPECIES", 100))

message(sprintf("DATA_DIR  = %s", Sys.getenv("DATA_DIR")))
message(sprintf("CACHE_DIR = %s", Sys.getenv("CACHE_DIR")))
message(sprintf("DEV_N_SPECIES = %s", Sys.getenv("DEV_N_SPECIES")))

# ---------------------------------------------------------------------------
# Step 1: Prepare data caches (reads from DATA_DIR, writes to CACHE_DIR)
# ---------------------------------------------------------------------------
message("\n=== Step 1: Prepare data ===")
source("_R/db_01_prepare-data.R")

# ---------------------------------------------------------------------------
# Step 2: Build KDE figure cache
# ---------------------------------------------------------------------------
message("\n=== Step 2: Build KDE cache ===")
source("_R/db_3_build-kde-cache.R")

# ---------------------------------------------------------------------------
# Step 3: Render (figure outputs land in data/ which _quarto.yml lists as
#          resources, so Quarto copies them to _site/data/ automatically)
# ---------------------------------------------------------------------------
message("\n=== Step 3: Render dashboard ===")

# Delete old-format cell_species JSON files so the render step regenerates them
# in the new compact format. Lookup files (var_codes, etc.) are always overwritten.
cs_dir    <- file.path("data", "cell_species")
old_jsons <- list.files(cs_dir, pattern = "^species_\\d{4}\\.json$", full.names = TRUE)
if (length(old_jsons)) {
  message(sprintf("Removing %d old cell_species JSON(s) for regeneration ...", length(old_jsons)))
  file.remove(old_jsons)
}

# Clear quarto freeze cache so all R chunks re-execute with fresh data
#unlink(".quarto/_freeze", recursive = TRUE)

exit_code <- system("quarto render .")
if (exit_code != 0) stop("quarto render failed (exit code ", exit_code, ")")

# ---------------------------------------------------------------------------
# Step 4: Publish to gh-pages
# ---------------------------------------------------------------------------
# Push _site/ as a clean orphan commit to origin/gh-pages.
# This bypasses quarto's branch management and works regardless of whether
# the remote branch exists or has stale large-file history.
message("\n=== Step 4: Publish to gh-pages ===")

# Get the remote URL from the main repo
remote_url <- trimws(system(
  "git -C ../.. remote get-url origin", intern = TRUE))
message("Remote: ", remote_url)

# Build a fresh single-commit repo inside _site/ and force-push
publish_cmds <- c(
  "cd _site",
  "git init -b gh-pages",
  "touch .nojekyll",
  "git add -A",
  'git -c user.email="build@local" -c user.name="build" commit -m "Deploy to GitHub Pages"',
  paste0('git push --force "', remote_url, '" gh-pages')
)
exit_code <- system(paste(publish_cmds, collapse = " && "))
if (exit_code != 0) stop("gh-pages push failed (exit code ", exit_code, ")")
if (exit_code != 0) stop("quarto publish failed (exit code ", exit_code, ")")

message("\n=== Done. Dashboard published to gh-pages. ===")
message("URL: https://cmerow.github.io/2025_Exposure/")
