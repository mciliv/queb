/**
 * Setup file for visual regression tests
 */

// Import the main setup first
require('./setup');

// Additional visual testing setup
const { toMatchImageSnapshot } = require('jest-image-snapshot');

expect.extend({ toMatchImageSnapshot });

// Configure image snapshot matching
beforeAll(() => {
  // Set up visual regression thresholds
  process.env.JEST_IMAGE_SNAPSHOT_TRACK_OBSOLETE = 'true';
});

// Mock canvas for visual tests
const { createCanvas } = require('canvas');
global.HTMLCanvasElement.prototype.getContext = function(type) {
  if (type === '2d') {
    return createCanvas(800, 600).getContext('2d');
  }
  return null;
};