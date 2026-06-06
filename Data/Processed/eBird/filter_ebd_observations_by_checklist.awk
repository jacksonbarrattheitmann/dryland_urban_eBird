
BEGIN {
  FS = "	"
  OFS = "	"
}

NR == FNR {
  if (FNR == 1 && $1 == "checklist_id") {
    next
  }
  ids[$1] = 1
  next
}

FNR == 1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col       = col["SAMPLING EVENT IDENTIFIER"]
  date_col     = col["OBSERVATION DATE"]
  common_col   = col["COMMON NAME"]
  sci_col      = col["SCIENTIFIC NAME"]
  tax_col      = col["TAXONOMIC ORDER"]
  category_col = col["CATEGORY"]
  count_col    = col["OBSERVATION COUNT"]
  lat_col      = col["LATITUDE"]
  lon_col      = col["LONGITUDE"]
  protocol_col = col["PROTOCOL NAME"]
  approved_col = col["APPROVED"]
  reviewed_col = col["REVIEWED"]

  missing = 0

  if (id_col == 0) {
    print "Missing SAMPLING EVENT IDENTIFIER" > "/dev/stderr"
    missing = 1
  }
  if (date_col == 0) {
    print "Missing OBSERVATION DATE" > "/dev/stderr"
    missing = 1
  }
  if (common_col == 0) {
    print "Missing COMMON NAME" > "/dev/stderr"
    missing = 1
  }
  if (sci_col == 0) {
    print "Missing SCIENTIFIC NAME" > "/dev/stderr"
    missing = 1
  }
  if (count_col == 0) {
    print "Missing OBSERVATION COUNT" > "/dev/stderr"
    missing = 1
  }

  if (missing == 1) {
    print "Columns found in EBD file:" > "/dev/stderr"
    for (i = 1; i <= NF; i++) {
      print i, $i > "/dev/stderr"
    }
    exit 1
  }

  print "checklist_id",
        "observation_date",
        "common_name",
        "scientific_name",
        "taxonomic_order",
        "category",
        "observation_count",
        "latitude",
        "longitude",
        "protocol_name",
        "approved",
        "reviewed"

  next
}

{
  checklist_id = $id_col

  if (checklist_id in ids) {
    print checklist_id,
          $date_col,
          $common_col,
          $sci_col,
          (tax_col > 0 ? $tax_col : ""),
          (category_col > 0 ? $category_col : ""),
          $count_col,
          (lat_col > 0 ? $lat_col : ""),
          (lon_col > 0 ? $lon_col : ""),
          (protocol_col > 0 ? $protocol_col : ""),
          (approved_col > 0 ? $approved_col : ""),
          (reviewed_col > 0 ? $reviewed_col : "")
  }
}

