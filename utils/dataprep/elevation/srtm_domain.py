def latitude_tiles(prefix, t_range):
    tiles = []
    for t in range(t_range[0], t_range[1] + 1):
        tile_part = "%s%02d" % (prefix, t)
        tiles.append(tile_part)
        
    return tiles  
  
def longitude_tiles(prefix, t_range):
    tiles = []
    for t in range(t_range[0], t_range[1] + 1):
        tile_part = "%s%03d" % (prefix, t)
        tiles.append(tile_part)
        
    return tiles
        
        
def all_tiles(lat_tiles, lon_tiles):
    tiles = []
    for lat in lat_tiles:
        for lon in lon_tiles:
            tiles.append("%s%s" % (lat, lon))
            
    return tiles

latitude_tiles_prefix = 'N'
latitude_tiles_range = (38, 42)

longitude_tiles_prefix = 'W'
longitude_tiles_range = (1,8)

tiles = all_tiles(latitude_tiles(latitude_tiles_prefix,latitude_tiles_range), 
    longitude_tiles(longitude_tiles_prefix, longitude_tiles_range))
