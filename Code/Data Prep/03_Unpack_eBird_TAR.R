# ============================================================
# 03_Unpack_eBird_TAR.R
#
# Purpose:
# Extract the eBird Sampling Event Data file from the downloaded
# July 2026 U.S. eBird TAR archive.
#
# The Sampling Event Data remains compressed as .txt.gz.
# Script 04 will stream this file directly without manually
# decompressing it.
#
# The full observation-level EBD remains inside the TAR because
# Script 05 streams it directly from the archive.
# ============================================================


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

ebird_raw_dir <- "Data/Raw/eBird"

ebird_tar <- file.path(
  ebird_raw_dir,
  "ebd_US_smp_relJul-2026.tar"
)

ebird_unzip_dir <- file.path(
  ebird_raw_dir,
  "Unpacked"
)

dir.create(
  ebird_unzip_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Expected Sampling Event Data member inside the TAR

sampling_member <- "ebd_US_relJul-2026_sampling.txt.gz"

sampling_out <- file.path(
  ebird_unzip_dir,
  sampling_member
)


# ------------------------------------------------------------
# 2. Confirm TAR archive exists
# ------------------------------------------------------------

if (!file.exists(ebird_tar)) {
  stop(
    "eBird TAR archive not found:\n",
    ebird_tar
  )
}


cat("\neBird TAR archive:\n")
print(ebird_tar)


cat("\nTAR archive size (GB):\n")

print(
  round(
    file.info(ebird_tar)$size / 1024^3,
    2
  )
)


# ------------------------------------------------------------
# 3. Inspect TAR contents
# ------------------------------------------------------------

tar_files <- utils::untar(
  tarfile = ebird_tar,
  list = TRUE
)


cat("\nFiles inside TAR archive:\n")
print(tar_files)


# Confirm Sampling Event Data is present

if (!sampling_member %in% tar_files) {
  stop(
    "Expected Sampling Event Data file was not found inside TAR:\n",
    sampling_member
  )
}


# Also confirm full EBD exists because Script 05 will need it

full_ebd_member <- "ebd_US_relJul-2026.txt.gz"


if (!full_ebd_member %in% tar_files) {
  stop(
    "Expected full EBD file was not found inside TAR:\n",
    full_ebd_member
  )
}


cat("\nRequired eBird files confirmed inside TAR:\n")

print(
  c(
    full_ebd_member,
    sampling_member
  )
)


# ------------------------------------------------------------
# 4. Remove any old Sampling Event Data extraction
# ------------------------------------------------------------

old_sampling_files <- list.files(
  ebird_unzip_dir,
  pattern = "sampling.*\\.txt(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)


if (length(old_sampling_files) > 0) {
  
  cat("\nRemoving old Sampling Event Data file(s):\n")
  print(old_sampling_files)
  
  file.remove(
    old_sampling_files
  )
}


# ------------------------------------------------------------
# 5. Extract Sampling Event Data only
# ------------------------------------------------------------

cat("\nExtracting Sampling Event Data...\n")


utils::untar(
  tarfile = ebird_tar,
  files = sampling_member,
  exdir = ebird_unzip_dir
)


# ------------------------------------------------------------
# 6. Confirm extraction succeeded
# ------------------------------------------------------------

if (!file.exists(sampling_out)) {
  stop(
    "Sampling Event Data extraction failed:\n",
    sampling_out
  )
}


sampling_size <- file.info(
  sampling_out
)$size


cat("\nExtracted Sampling Event Data:\n")
print(sampling_out)


cat("\nSampling Event Data size (GB):\n")

print(
  round(
    sampling_size / 1024^3,
    3
  )
)


if (is.na(sampling_size) || sampling_size == 0) {
  stop(
    "Sampling Event Data file is empty. Do not continue to Script 04."
  )
}


# ------------------------------------------------------------
# 7. Confirm Sampling Event Data header
# ------------------------------------------------------------

library(data.table)


sampling_header <- fread(
  cmd = paste(
    "gzip -cd",
    shQuote(sampling_out)
  ),
  sep = "\t",
  quote = "",
  nrows = 0,
  showProgress = FALSE
)


cat("\nSampling Event Data columns:\n")
print(names(sampling_header))


required_columns <- c(
  "SAMPLING EVENT IDENTIFIER",
  "OBSERVATION DATE",
  "LATITUDE",
  "LONGITUDE",
  "ALL SPECIES REPORTED",
  "DURATION MINUTES",
  "EFFORT DISTANCE KM",
  "NUMBER OBSERVERS"
)


missing_columns <- setdiff(
  required_columns,
  names(sampling_header)
)


if (length(missing_columns) > 0) {
  
  stop(
    "Sampling Event Data is missing required columns:\n",
    paste(
      missing_columns,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 8. Final confirmation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 03 COMPLETE\n")
cat("============================================================\n\n")

cat("Sampling Event Data successfully extracted and verified:\n")
cat(sampling_out, "\n")

cat("\nFull observation EBD remains inside TAR for Script 05:\n")
cat(full_ebd_member, "\n")

cat("\nReady for Script 04.\n")