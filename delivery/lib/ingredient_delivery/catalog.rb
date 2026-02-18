module IngredientDelivery
  class Catalog
    attr_reader :ingredients

    def initialize(data_path:)
      @ingredients = JSON.parse(File.read(data_path), symbolize_names: true)
    end
  end
end
