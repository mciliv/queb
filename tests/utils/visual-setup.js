/**
 * Visual test setup
 * Additional setup for visual regression tests
 */

const fs = require('fs');
const path = require('path');

// Directory for visual test snapshots
const SNAPSHOTS_DIR = path.join(__dirname, '../__snapshots__/visual');
const DIFF_DIR = path.join(__dirname, '../__snapshots__/diff');

beforeAll(() => {
  // Ensure snapshot directories exist
  [SNAPSHOTS_DIR, DIFF_DIR].forEach(dir => {
    if (!fs.existsSync(dir)) {
      fs.mkdirSync(dir, { recursive: true });
    }
  });
  
  // Set up visual regression configuration
  global.VISUAL_CONFIG = {
    snapshotsDir: SNAPSHOTS_DIR,
    diffDir: DIFF_DIR,
    threshold: parseFloat(process.env.VISUAL_REGRESSION_THRESHOLD || '0.01'),
    updateSnapshots: process.env.UPDATE_SNAPSHOTS === 'true'
  };
});

// Helper to compare images
global.compareImages = async (actualImage, expectedImage, diffImage) => {
  // This would use a library like pixelmatch or jest-image-snapshot
  // For now, we'll create a simple mock
  const { PNG } = require('pngjs');
  const pixelmatch = require('pixelmatch');
  
  const img1 = PNG.sync.read(fs.readFileSync(actualImage));
  const img2 = PNG.sync.read(fs.readFileSync(expectedImage));
  
  const { width, height } = img1;
  const diff = new PNG({ width, height });
  
  const numDiffPixels = pixelmatch(
    img1.data,
    img2.data,
    diff.data,
    width,
    height,
    { threshold: global.VISUAL_CONFIG.threshold }
  );
  
  if (numDiffPixels > 0) {
    fs.writeFileSync(diffImage, PNG.sync.write(diff));
  }
  
  return {
    diffPixels: numDiffPixels,
    diffPercentage: (numDiffPixels / (width * height)) * 100
  };
};

// Helper to capture screenshot with consistent settings
global.captureScreenshot = async (page, name, options = {}) => {
  const screenshotPath = path.join(SNAPSHOTS_DIR, `${name}.png`);
  
  // Wait for animations to complete
  await page.waitForTimeout(options.waitForAnimations || 500);
  
  // Hide elements that might cause flakiness (e.g., cursors, timestamps)
  if (options.hideSelectors) {
    for (const selector of options.hideSelectors) {
      await page.addStyleTag({
        content: `${selector} { visibility: hidden !important; }`
      });
    }
  }
  
  // Take screenshot
  await page.screenshot({
    path: screenshotPath,
    fullPage: options.fullPage !== false,
    clip: options.clip
  });
  
  return screenshotPath;
};

// Helper for visual regression testing
global.expectVisualMatch = async (actualPath, testName) => {
  const expectedPath = path.join(SNAPSHOTS_DIR, 'baseline', `${testName}.png`);
  const diffPath = path.join(DIFF_DIR, `${testName}-diff.png`);
  
  // If updating snapshots or baseline doesn't exist, save as new baseline
  if (global.VISUAL_CONFIG.updateSnapshots || !fs.existsSync(expectedPath)) {
    const baselineDir = path.dirname(expectedPath);
    if (!fs.existsSync(baselineDir)) {
      fs.mkdirSync(baselineDir, { recursive: true });
    }
    fs.copyFileSync(actualPath, expectedPath);
    console.log(`Saved new baseline snapshot: ${testName}`);
    return;
  }
  
  // Compare images
  const result = await global.compareImages(actualPath, expectedPath, diffPath);
  
  if (result.diffPercentage > 0.1) { // 0.1% threshold
    throw new Error(
      `Visual regression detected for ${testName}: ` +
      `${result.diffPercentage.toFixed(2)}% difference ` +
      `(${result.diffPixels} pixels). ` +
      `Diff saved to: ${diffPath}`
    );
  }
};

// Helper to test responsive design
global.testResponsive = async (page, url, testName) => {
  const viewports = [
    { name: 'mobile', width: 375, height: 667 },
    { name: 'tablet', width: 768, height: 1024 },
    { name: 'desktop', width: 1920, height: 1080 }
  ];
  
  for (const viewport of viewports) {
    await page.setViewport(viewport);
    await page.goto(url, { waitUntil: 'networkidle0' });
    
    const screenshotPath = await global.captureScreenshot(
      page,
      `${testName}-${viewport.name}`
    );
    
    await global.expectVisualMatch(
      screenshotPath,
      `${testName}-${viewport.name}`
    );
  }
};

// Custom matchers for visual tests
expect.extend({
  async toMatchVisualSnapshot(page, name, options = {}) {
    try {
      const screenshotPath = await global.captureScreenshot(page, name, options);
      await global.expectVisualMatch(screenshotPath, name);
      
      return {
        pass: true,
        message: () => `Visual snapshot matched for ${name}`
      };
    } catch (error) {
      return {
        pass: false,
        message: () => error.message
      };
    }
  }
});
