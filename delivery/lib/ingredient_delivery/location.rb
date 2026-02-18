module IngredientDelivery
  module Location
    EARTH_RADIUS_KM = 6371.0

    def self.distance_km(lat1, lng1, lat2, lng2)
      lat1_rad = degrees_to_radians(lat1)
      lat2_rad = degrees_to_radians(lat2)
      delta_lat = degrees_to_radians(lat2 - lat1)
      delta_lng = degrees_to_radians(lng2 - lng1)

      a = Math.sin(delta_lat / 2)**2 +
          Math.cos(lat1_rad) * Math.cos(lat2_rad) *
          Math.sin(delta_lng / 2)**2
      c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a))

      (EARTH_RADIUS_KM * c).round(2)
    end

    def self.degrees_to_radians(degrees)
      degrees * Math::PI / 180.0
    end
  end
end
