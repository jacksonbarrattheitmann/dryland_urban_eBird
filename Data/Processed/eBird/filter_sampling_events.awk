
BEGIN {
  FS = "	"
  OFS = ","
}

NR == 1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col       = col["SAMPLING EVENT IDENTIFIER"]
  date_col     = col["OBSERVATION DATE"]
  lat_col      = col["LATITUDE"]
  lon_col      = col["LONGITUDE"]
  protocol_col = col["PROTOCOL TYPE"]
  complete_col = col["ALL SPECIES REPORTED"]
  duration_col = col["DURATION MINUTES"]
  distance_col = col["EFFORT DISTANCE KM"]
  observer_col = col["NUMBER OBSERVERS"]

  if (
    id_col == 0 ||
    date_col == 0 ||
    lat_col == 0 ||
    lon_col == 0 ||
    protocol_col == 0 ||
    complete_col == 0 ||
    duration_col == 0 ||
    distance_col == 0 ||
    observer_col == 0
  ) {
    print "ERROR: One or more required columns were not found." > "/dev/stderr"
    print "Columns found in file:" > "/dev/stderr"
    for (i = 1; i <= NF; i++) {
      print i, $i > "/dev/stderr"
    }
    exit 1
  }

  print "checklist_id","observation_date","latitude","longitude","protocol_type","all_species_reported","duration_minutes","effort_distance_km","number_observers"
  next
}

{
  checklist_id = $id_col
  obs_date     = $date_col
  lat_raw      = $lat_col
  lon_raw      = $lon_col
  lat          = lat_raw + 0
  lon          = lon_raw + 0
  protocol     = $protocol_col
  complete     = $complete_col
  duration     = $duration_col + 0
  distance     = $distance_col
  observers    = $observer_col + 0

  if (
    obs_date >= "2015-01-01" &&
    obs_date <= "2025-12-31" &&
    protocol == "Stationary" &&
    (complete == "1" || complete == "TRUE" || complete == "true" || complete == "T") &&
    observers == 1 &&
    duration >= 5 &&
    duration <= 240 &&
    lat_raw != "" &&
    lon_raw != "" &&
    lat >= ymin &&
    lat <= ymax &&
    lon >= xmin &&
    lon <= xmax
  ) {
    print checklist_id,obs_date,lat,lon,protocol,complete,duration,distance,observers
  }
}

