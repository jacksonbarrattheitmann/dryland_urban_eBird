# ============================================================
# 03_Unpack_eBird_TAR.R
#
# Purpose:
# Unpack the downloaded eBird .tar file without opening it
# manually on the desktop.
# I need to do this because its so big I cannot open it on its own.
# ============================================================

library(fs)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

ebird_raw_dir <- "Data/Raw/eBird"
ebird_tar <- file.path(ebird_raw_dir, "ebd_US_smp_relApr-2026.tar")

ebird_unzip_dir <- file.path(ebird_raw_dir, "Unpacked")
dir_create(ebird_unzip_dir)

# ------------------------------------------------------------
# 2. Inspect files inside the TAR archive
# ------------------------------------------------------------

tar_files <- utils::untar(
  tarfile = ebird_tar,
  list = TRUE
)

print(tar_files)

# ------------------------------------------------------------
# 3. Extract TAR archive
# ------------------------------------------------------------

utils::untar(
  tarfile = ebird_tar,
  exdir = ebird_unzip_dir
)

# ------------------------------------------------------------
# 4. List extracted files
# ------------------------------------------------------------

extracted_files <- dir(
  ebird_unzip_dir,
  recursive = TRUE,
  full.names = TRUE
)

print(extracted_files)