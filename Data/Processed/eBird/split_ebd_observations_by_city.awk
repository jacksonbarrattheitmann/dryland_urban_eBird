
BEGIN {
  FS = "	"
  OFS = "	"
}

NR == FNR {
  if (FNR == 1 && $1 == "checklist_id") {
    next
  }

  checklist_to_city_file[$1] = $2
  next
}

FNR == 1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col = col["SAMPLING EVENT IDENTIFIER"]
  date_col = col["OBSERVATION DATE"]
  locality_id_col = col["LOCALITY ID"]
  locality_col = col["LOCALITY"]
  common_col = col["COMMON NAME"]
  sci_col = col["SCIENTIFIC NAME"]
  tax_col = col["TAXONOMIC ORDER"]
  category_col = col["CATEGORY"]
  count_col = col["OBSERVATION COUNT"]
  lat_col = col["LATITUDE"]
  lon_col = col["LONGITUDE"]
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

  if (locality_id_col == 0) {
    print "Missing LOCALITY ID" > "/dev/stderr"
    missing = 1
  }

  if (locality_col == 0) {
    print "Missing LOCALITY" > "/dev/stderr"
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

  header = "checklist_id"
  header = header OFS "observation_date"
  header = header OFS "locality_id"
  header = header OFS "locality"
  header = header OFS "common_name"
  header = header OFS "scientific_name"
  header = header OFS "taxonomic_order"
  header = header OFS "category"
  header = header OFS "observation_count"
  header = header OFS "latitude"
  header = header OFS "longitude"
  header = header OFS "protocol_name"
  header = header OFS "approved"
  header = header OFS "reviewed"

  next
}

{
  checklist_id = $id_col

  if (checklist_id in checklist_to_city_file) {
    outfile = outdir "/" checklist_to_city_file[checklist_id] ".tsv"

    if (!(outfile in wrote_header)) {
      print header > outfile
      wrote_header[outfile] = 1
      close(outfile)
    }

    locality_id_val = ""
    locality_val = ""
    tax_val = ""
    category_val = ""
    lat_val = ""
    lon_val = ""
    protocol_val = ""
    approved_val = ""
    reviewed_val = ""

    if (locality_id_col > 0) {
      locality_id_val = $locality_id_col
    }

    if (locality_col > 0) {
      locality_val = $locality_col
    }

    if (tax_col > 0) {
      tax_val = $tax_col
    }

    if (category_col > 0) {
      category_val = $category_col
    }

    if (lat_col > 0) {
      lat_val = $lat_col
    }

    if (lon_col > 0) {
      lon_val = $lon_col
    }

    if (protocol_col > 0) {
      protocol_val = $protocol_col
    }

    if (approved_col > 0) {
      approved_val = $approved_col
    }

    if (reviewed_col > 0) {
      reviewed_val = $reviewed_col
    }

    row = checklist_id
    row = row OFS $date_col
    row = row OFS locality_id_val
    row = row OFS locality_val
    row = row OFS $common_col
    row = row OFS $sci_col
    row = row OFS tax_val
    row = row OFS category_val
    row = row OFS $count_col
    row = row OFS lat_val
    row = row OFS lon_val
    row = row OFS protocol_val
    row = row OFS approved_val
    row = row OFS reviewed_val

    print row >> outfile
    close(outfile)
  }
}

