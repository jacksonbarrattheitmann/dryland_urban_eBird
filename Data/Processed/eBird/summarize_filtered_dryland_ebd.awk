
BEGIN {
  FS = "\t"
  OFS = "\t"
}

FNR == 1 {

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  city_col = col["city_file_stub"]
  id_col = col["checklist_id"]
  date_col = col["observation_date"]
  sci_col = col["scientific_name"]
  category_col = col["category"]
  protocol_col = col["protocol_name"]

  # Keep the full condition on one line for macOS / BSD awk.
  # BSD awk can reject a newline immediately after `if (`.
  if (city_col == 0 || id_col == 0 || date_col == 0 || sci_col == 0 || category_col == 0 || protocol_col == 0) {
    print "Filtered EBD is missing required summary columns." > "/dev/stderr"
    exit 1
  }

  next
}

{

  checklist_id = $id_col
  protocol_name = $protocol_col

  if (protocol_name != "") {
    print checklist_id, protocol_name >> protocol_pairs
  }

  if ($category_col == "species" && $sci_col != "") {

    scientific_name = $sci_col
    city_stub = $city_col
    year = substr($date_col, 1, 4)

    print checklist_id, scientific_name >> checklist_pairs
    print city_stub, year, scientific_name >> city_year_pairs
  }
}

