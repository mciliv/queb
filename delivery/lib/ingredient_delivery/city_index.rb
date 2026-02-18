module IngredientDelivery
  class CityIndex
    def initialize(data_path:)
      @cities = JSON.parse(File.read(data_path), symbolize_names: true)
    end

    def fetch(city_name)
      key = normalize(city_name)
      city = @cities.find { |entry| normalize(entry[:name]) == key }
      raise ArgumentError, "Unknown city: #{city_name}" unless city

      { name: city[:name], lat: city[:lat], lng: city[:lng] }
    end

    private

    def normalize(value)
      value.to_s.strip.downcase
    end
  end
end
