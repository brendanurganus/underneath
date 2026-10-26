# Underneath

An interactive 3D explorer of earthquakes and the inside of the Earth.

- Every magnitude 4.5+ earthquake of the past year from the USGS catalog, plus every M8.5+ quake since 1900
- Slice the planet open to see how deep quakes happen
- **Wave Race:** watch P, S and surface waves cross a cutaway Earth, alongside real seismograms from global stations
- **Year in 30 seconds:** replay a year of quakes and watch the energy add up
- Tectonic plate boundaries and speeds, and 1,214 volcanoes with current US alert levels

Not a warning or alert system.

## Run it

It's a static site. Serve the folder with any web server (opening `index.html` as a file won't load the data):

```
python3 -m http.server
```

Then open http://localhost:8000.

## Data sources

- Earthquakes: USGS ANSS Comprehensive Catalog (ComCat)
- Seismograms: EarthScope (IRIS) FDSN data services
- Plate boundaries: Bird (2003), PB2002
- Volcanoes: Smithsonian Global Volcanism Program; alert levels from the USGS Volcano Hazards Program
- Relief: NOAA ETOPO1
- Earth model and travel times: iasp91 via ObsPy TauP

## Updating the data

`scripts/refresh.py` rebuilds `data/` from the live sources:

```
python3 scripts/refresh.py . out
```

The other scripts in `scripts/` rebuild the static files (relief texture, plates, travel-time tables) from scratch.
