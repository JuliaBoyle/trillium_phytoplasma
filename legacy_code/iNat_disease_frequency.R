library(tidyverse)
library(sf)
library(mapview)
library(viridis)
library(rnaturalearth)

# Set seed for subsampling by user login
set.seed(1)

dat = read_csv("iNat_alltrilliums_2023_07_05.csv") 


# Check out distribution of observatiosn across users
dat %>% 
  group_by(user_login) %>% 
  summarize(n = n()) %>% 
  arrange(-n)

dat %>% distinct(user_login)
# 16158 distinct observers

dat %>% 
  group_by(user_login) %>% 
  summarize(n = n()) %>% 
  filter(n>15)
# 194 users had >15 observations

dat %>% 
  group_by(user_login) %>% 
  summarize(n_obs = n()) %>% 
  group_by(n_obs) %>% 
  summarise(n_users = n()) %>% 
  ggplot() +
    geom_point(aes(x = n_obs, y = n_users))
 

dat_cleaned = dat %>% 
  # Adding positional accuracy filter to see if that helps deal with observations in lakes
  filter(!is.na(latitude), !is.na(longitude), !is.na(public_positional_accuracy)) %>% 
  filter(public_positional_accuracy < 10000) %>% 
  select(latitude, longitude, infected = `field:phytoplasma infection?`, observed_on, user_login, id) %>% 
  mutate(infected = if_else(is.na(infected), "No", infected)) %>% 
  # Drop european and puerto rico observations
  filter(latitude > 20, longitude < -65) %>% 
  # Group by user and subset to 15 observations per user
  group_by(user_login) %>% 
  slice_sample(n = 15, replace = FALSE)


# inat<-dat %>% mutate(country_code=ifelse(grepl("Canada", place_country_name), "CAN", "USA"))
# inat_clean <- clean_coordinates(x = dat, 
#                                 lon = "longitude", 
#                                 lat = "latitude",
#                                 countries = "country_code",
#                                 species = "scientific_name",
#                                 tests = c("centroids",
#                                           "equal", "zeros", "institutions", "duplicates")) 
# inat_clean=inat_clean %>% filter(.summary==TRUE) #only include observations that have no suspect issues

# Quick look at infections over time
ggplot(dat_cleaned) +
  geom_histogram(aes(x = observed_on, fill = infected))

# Convert to simple features object
dat_sf = st_as_sf(dat_cleaned, coords = c("longitude", "latitude")) %>% 
  st_set_crs("EPSG:4326")

# Quick look
plot(dat_sf)

# Reproject to equal area projection
dat_ea = st_transform(dat_sf, crs = "ESRI:102003")

# Quick look
plot(dat_ea)

# Interactive map, can check for points in lakes etc. 
mapview(dat_ea)

# Make a grid of polygons that encompass the observation points
# Choosing an arbitrary 5 km grid size
area_fishnet_grid = st_make_grid(dat_ea, c(5e4, 5e4), what = "polygons", square = FALSE)

# Convert to sf and add grid ID
fishnet_grid_sf = st_sf(area_fishnet_grid) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(area_fishnet_grid)))

# Quick look at the grid
mapview(fishnet_grid_sf)

# Assign observations to grid cells
grid_assignments = st_intersection(fishnet_grid_sf, dat_ea)

# Calculate frequencies per grid cell
freq_per_poly = grid_assignments %>% 
  group_by(grid_id) %>% 
  summarize(n = n(), n_infected = sum(infected == "Yes")) %>% 
  mutate(freq_infected = n_infected/n)

# Rejoin these calculated frequencies to polygons for plotting
freqs_with_polys = st_join(fishnet_grid_sf, freq_per_poly) %>% 
  filter(!is.na(n))

# Quick look 
plot(freqs_with_polys)

world = ne_countries(scale = "large", returnclass = "sf", continent = "north america")
lakes = ne_download(scale = "medium", type = "lakes", category = "physical")

ggplot() +
  geom_sf(data = world, fill = "white") +
  geom_sf(data = lakes, color = "grey", fill = "grey") +
  # geom_sf(data = dat_sf, alpha = 01, size = 0.1, color = "grey20") +
  geom_sf(data = freqs_with_polys, aes(fill = freq_infected), alpha = 0.8, color = "lightgrey") +
  # scale_fill_viridis(direction = -1, option = "G") +
  scale_fill_gradientn(colours = c("white", "lightblue", "navy"), values = c(0,0.000001,1)) +
  coord_sf(xlim = c(-102.15, -68), ylim = c(30, 55), expand = FALSE) +
  theme_bw()
      
ggplot() +
  geom_sf(data = world, fill = "white") +
  geom_sf(data = lakes, color = "grey", fill = "grey") +
  # geom_sf(data = dat_sf, alpha = 01, size = 0.1, color = "grey20") +
  geom_sf(data = freqs_with_polys, aes(fill = n), alpha = 0.8, color = "lightgrey") +
  # scale_fill_viridis(direction = -1, option = "G") +
  scale_fill_gradientn(colours = c("lightblue", "navy"), values = c(0,1)) +
  coord_sf(xlim = c(-102.15, -68), ylim = c(30, 55), expand = FALSE) +
  theme_bw()

ggplot(freqs_with_polys, aes(x = n, y = freq_infected)) +
  geom_point() +
  geom_smooth()
  

