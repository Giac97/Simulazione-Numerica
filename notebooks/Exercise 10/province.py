from geopy.geocoders import Nominatim

# Lista di capoluoghi di provincia italiani
capoluoghi = []

with open("province.txt", 'r') as file:
    capoluoghi = file.read().splitlines()

file.close()
# Inizializza il geolocalizzatore
geolocator = Nominatim(user_agent="myGeocoder")

# Crea una lista di tuple con il nome del capoluogo e le sue coordinate
coordinate_capoluoghi = [(capoluogo, geolocator.geocode(capoluogo + ", Italy")) for capoluogo in capoluoghi]

# Stampa le coordinate
with open("capoluoghi.txt", 'w') as file:
    for capoluogo, location in coordinate_capoluoghi:
        if location:
            file.write(f"{location.longitude} {location.latitude}\n")

file.close()
