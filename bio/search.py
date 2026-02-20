import os
import requests
import json

# Your API key from USDA FoodData Central
api_key = os.environ.get('FDC_API_KEY')
if not api_key:
    raise ValueError("API Key not found")

# Base URL for the FoodData Central API
base_url = "https://api.nal.usda.gov/fdc/v1/foods/search"

# Function to search for a food item
def search_food(food_name):
    params = {
        "api_key": api_key,
        "query": food_name,
        "dataType": ["Survey (FNDDS)", "Foundation", "SR Legacy"],  # Types of food data to search for
        "pageSize": 1  # We're only interested in the first result for simplicity
    }
    
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = response.json()
        if data.get('foods'):
            return data['foods'][0]
    return None

# List of food items
food_items = [
    "Dried shiitake", "Dried maitake", "Avocado", "Olive", "Broccoli", "Broccoli Rabe",
    "Cauliflower", "Spinach", "Arugula", "Watercress", "Kale", "Fennel", "Purslane", 
    "Chard", "Celery", "Red cabbage", "Asparagus", "Carrot", "Beet", "Bittermelon",
    "Tomatoes", "Small Tomatoes", "Heirloom Tomatoes", "Sweet potato", "Japanese Sweet potato",
    "Potato", "Lime", "Shallot", "Chile", "Onion", "Green Onion", "Red Onion", "Yellow Onion",
    "Garlic", "Caper", "Pickled Caper", "Non-pareil Pickled Caper", "Turmeric", "Ginger", 
    "Ginseng", "Cilantro", "Parsley", "Pomegranate", "Kiwi", "Apple", "Eggs", "Caviar"
]

# Dictionary to hold nutritional data for all items
nutritional_data = {}

# Fetch data for each item
for item in food_items:
    food_info = search_food(item)
    if food_info:
        # Here we're extracting some common nutritional values but you can customize this based on what's available
        nutritional_data[item] = {
            "description": food_info.get('description'),
            "energy": next((nutrient for nutrient in food_info.get('foodNutrients', []) if nutrient.get('nutrientName') == 'Energy'), None)['value'] if food_info.get('foodNutrients') else None,
            "protein": next((nutrient for nutrient in food_info.get('foodNutrients', []) if nutrient.get('nutrientName') == 'Protein'), None)['value'] if food_info.get('foodNutrients') else None,
            "fat": next((nutrient for nutrient in food_info.get('foodNutrients', []) if nutrient.get('nutrientName') == 'Total lipid (fat)'), None)['value'] if food_info.get('foodNutrients') else None,
            "carbohydrates": next((nutrient for nutrient in food_info.get('foodNutrients', []) if nutrient.get('nutrientName') == 'Carbohydrate, by difference'), None)['value'] if food_info.get('foodNutrients') else None
        }
    else:
        nutritional_data[item] = {"error": "No data found for this item."}

# Save the data to a JSON file
with open('usda_nutritional_info.json', 'w') as json_file:
    json.dump(nutritional_data, json_file, indent=4)

print("Nutritional data has been saved to 'usda_nutritional_info.json'")