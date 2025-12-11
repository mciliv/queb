// Minimal screenshot service stub used for tests. The real implementation
// is mocked in unit and integration tests; this file exists to satisfy module resolution.

class ScreenshotService {
  async captureApp() {
    return { success: true, filename: 'stub.png' };
  }
  async captureWithInput() {
    return { success: true, filename: 'stub.png' };
  }
  async captureAnalysis() {
    return { success: true, filename: 'stub.png' };
  }
  async listScreenshots() {
    return [];
  }
  async getScreenshotPath() {
    return '/tmp/stub.png';
  }
  async cleanupOldScreenshots() {
    return { cleaned: 0 };
  }
}

module.exports = ScreenshotService;
