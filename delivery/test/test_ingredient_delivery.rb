require "minitest/autorun"
require_relative "../lib/ingredient_delivery"

class IngredientDeliveryTest < Minitest::Test
  def setup
    data_root = File.expand_path("../data", __dir__)
    @catalog = IngredientDelivery.load_catalog(data_path: File.join(data_root, "ingredients.json"))
    @city_index = IngredientDelivery.load_city_index(data_path: File.join(data_root, "cities.json"))
    @matcher = IngredientDelivery::Matcher.new(catalog: @catalog, city_index: @city_index)
  end

  def test_best_for_city_returns_limited_results
    results = @matcher.best_for_city("San Francisco", limit: 3)

    assert_equal 3, results.length
    assert results.all? { |entry| entry[:user_location][:name] == "San Francisco" }
  end

  def test_results_are_sorted_by_score_desc
    results = @matcher.best_for_city("Seattle", limit: 5)
    scores = results.map { |entry| entry[:score] }

    assert_equal scores.sort.reverse, scores
  end

  def test_best_for_location_includes_delivery_focus
    results = @matcher.best_for_location(lat: 40.7128, lng: -74.006, limit: 2)

    assert results.all? { |entry| entry[:delivery_focus] == "highest_quality" }
  end
end
