
BEGIN {
  FS = "	"
  OFS = "	"
}

FNR == 1 {

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col = col["SAMPLING EVENT IDENTIFIER"]
  date_col = col["OBSERVATION DATE"]
  lat_col = col["LATITUDE"]
  lon_col = col["LONGITUDE"]
  complete_col = col["ALL SPECIES REPORTED"]
  duration_col = col["DURATION MINUTES"]
  distance_col = col["EFFORT DISTANCE KM"]
  observer_col = col["NUMBER OBSERVERS"]

  locality_id_col = col["LOCALITY ID"]
  locality_col = col["LOCALITY"]

  print "checklist_id", "observation_date", "latitude", "longitude", "all_species_reported", "duration_minutes", "effort_distance_km", "number_observers", "locality_id", "locality"

  next
}

{

  date_val = $date_col
  lat_val = $lat_col
  lon_val = $lon_col
  duration_val = $duration_col
  observer_val = $observer_col
  complete_val = $complete_col

  if (date_val < start_date || date_val > end_date) next
  if (complete_val != 1) next
  if (observer_val != 1) next
  if (duration_val < 5 || duration_val > 240) next

  if (lat_val == "") next
  if (lon_val == "") next

  if ((lat_val + 0) < ymin || (lat_val + 0) > ymax) next
  if ((lon_val + 0) < xmin || (lon_val + 0) > xmax) next

  locality_id_val = ""
  locality_val = ""
  distance_val = ""

  if (locality_id_col > 0) {
    locality_id_val = $locality_id_col
  }

  if (locality_col > 0) {
    locality_val = $locality_col
  }

  if (distance_col > 0) {
    distance_val = $distance_col
  }

  print $id_col, date_val, lat_val, lon_val, complete_val, duration_val, distance_val, observer_val, locality_id_val, locality_val
}

