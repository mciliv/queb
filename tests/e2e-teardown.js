/**
 * Global teardown for E2E tests
 */

module.exports = async () => {
  console.log('Tearing down E2E test environment...');

  // Clean up test database or services
  // if (global.testServer) {
  //   await global.testServer.close();
  // }

  // Clean up any test data
  // await cleanupTestData();
};