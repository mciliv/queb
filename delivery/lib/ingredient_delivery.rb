require "json"
require_relative "ingredient_delivery/catalog"
require_relative "ingredient_delivery/city_index"
require_relative "ingredient_delivery/location"
require_relative "ingredient_delivery/matcher"

module IngredientDelivery
  def self.load_catalog(data_path:)
    Catalog.new(data_path: data_path)
  end

  def self.load_city_index(data_path:)
    CityIndex.new(data_path: data_path)
  end
end
