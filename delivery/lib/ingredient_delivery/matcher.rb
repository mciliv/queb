module IngredientDelivery
  class Matcher
    DEFAULT_LIMIT = 5

    def initialize(catalog:, city_index:)
      @catalog = catalog
      @city_index = city_index
    end

    def best_for_city(city_name, limit: DEFAULT_LIMIT)
      city = @city_index.fetch(city_name)
      best_for_location(lat: city[:lat], lng: city[:lng], limit: limit, location_name: city[:name])
    end

    def best_for_location(lat:, lng:, limit: DEFAULT_LIMIT, location_name: nil)
      # Block: best_for_location.
      # Refined prompt: build a Ruby app that selects the highest-quality ingredients
      # for a user's location using clear, explainable scoring and proximity-aware
      # ranking. Prompt: "in this project, create a ruby app to deliver highest quality
      # ingredients based on user's loc".
      scored = @catalog.ingredients.map do |ingredient|
        distance = Location.distance_km(lat, lng, ingredient[:lat], ingredient[:lng])
        score = score_ingredient(ingredient, distance)
        ingredient.merge(score: score, distance_km: distance)
      end

      sorted = scored.sort_by { |entry| -entry[:score] }
      limited = sorted.first(limit)
      add_context(limited, lat: lat, lng: lng, location_name: location_name)
    end

    private

    def score_ingredient(ingredient, distance_km)
      quality = ingredient[:quality_score].to_f
      supplier = ingredient[:supplier_rating].to_f freshness_score = [0.0, 100.0 - ingredient[:freshness_hours].to_f * 2.0].max
      proximity_score = [0.0, 100.0 - distance_km * 2.0].max

      (quality * 0.55 + supplier * 0.2 + freshness_score * 0.15 + proximity_score * 0.1).round(2)
    end

    def add_context(results, lat:, lng:, location_name:)
      results.map do |entry|
        entry.merge(
          delivery_focus: "highest_quality",
          user_location: {
            name: location_name,
            lat: lat,
            lng: lng
          }
        )
      end
    end
  end
end
