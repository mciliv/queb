const https = require('https');
const http = require('http');
const fs = require('fs');
const path = require('path');

/**
 * Generic test image downloader
 * Usage: node download-test-image.js <name> <image_url>
 * Example: node download-test-image.js egg "https://images.pexels.com/photos/1556707/pexels-photo-1556707.jpeg?auto=compress&cs=tinysrgb&w=800"
 */

function downloadTestImage(name, imageUrl) {
  // Ensure images directory exists
  const imagesDir = path.join(__dirname, 'public', 'images');
  if (!fs.existsSync(imagesDir)) {
    fs.mkdirSync(imagesDir, { recursive: true });
  }

  // Determine file extension
  let ext = 'jpg';
  const urlParts = imageUrl.split('.');
  if (urlParts.length > 1) {
    const possibleExt = urlParts[urlParts.length - 1].split('?')[0].toLowerCase();
    if (['jpg', 'jpeg', 'png', 'webp', 'gif'].includes(possibleExt)) {
      ext = possibleExt;
    }
  }

  const outputPath = path.join(imagesDir, `${name}.${ext}`);
  const file = fs.createWriteStream(outputPath);

  // Choose protocol based on URL
  const protocol = imageUrl.startsWith('https') ? https : http;

  // Download the image
  protocol.get(imageUrl, (response) => {
    response.pipe(file);
    
    file.on('finish', () => {
      file.close();
      console.log(`✅ Successfully downloaded ${name} image`);
      console.log(`📸 Saved to: ${outputPath}`);
      console.log(`📊 Size: ${fs.statSync(outputPath).size.toLocaleString()} bytes`);
    });
  }).on('error', (err) => {
    fs.unlink(outputPath, () => {});
    console.error('❌ Error downloading image:', err.message);
    process.exit(1);
  });
}

// Parse command line arguments
const args = process.argv.slice(2);

if (args.length === 0) {
  // Default to egg image
  const name = 'egg';
  const url = 'https://images.pexels.com/photos/1556707/pexels-photo-1556707.jpeg?auto=compress&cs=tinysrgb&w=800';
  downloadTestImage(name, url);
} else if (args.length === 2) {
  const [name, url] = args;
  downloadTestImage(name, url);
} else {
  console.log('Usage: node download-test-image.js <name> <image_url>');
  console.log('Example: node download-test-image.js egg "https://example.com/egg.jpg"');
  console.log('\nOr run without arguments to download default egg image');
  process.exit(1);
}
