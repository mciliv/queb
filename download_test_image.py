import requests
import os
import sys

# Generic test image downloader
# Usage: python download_test_image.py <name> <image_url>
# Example: python download_test_image.py egg "https://images.pexels.com/photos/1556707/pexels-photo-1556707.jpeg?auto=compress&cs=tinysrgb&w=800"

def download_test_image(name, image_url):
    """Download a test image for testing"""
    # Create images directory if it doesn't exist
    os.makedirs("public/images", exist_ok=True)
    
    # Determine file extension from URL or default to jpg
    ext = 'jpg'
    if '.' in image_url:
        url_ext = image_url.split('.')[-1].split('?')[0].lower()
        if url_ext in ['jpg', 'jpeg', 'png', 'webp', 'gif']:
            ext = url_ext
    
    output_path = f"public/images/{name}.{ext}"
    
    # Download the image
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        response = requests.get(image_url, headers=headers)
        
        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
            print(f"✅ Successfully downloaded {name} image")
            print(f"📸 Image size: {len(response.content):,} bytes")
            print(f"📍 Saved to: {output_path}")
            return True
        else:
            print(f"❌ Failed to download image: {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ Error downloading image: {e}")
        return False

if __name__ == "__main__":
    # Default to egg image if no args provided
    if len(sys.argv) == 1:
        name = "egg"
        url = "https://images.pexels.com/photos/1556707/pexels-photo-1556707.jpeg?auto=compress&cs=tinysrgb&w=800"
    elif len(sys.argv) == 3:
        name = sys.argv[1]
        url = sys.argv[2]
    else:
        print("Usage: python download_test_image.py <name> <image_url>")
        print("Example: python download_test_image.py egg 'https://example.com/egg.jpg'")
        sys.exit(1)
    
    download_test_image(name, url)
