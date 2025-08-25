// Google Cloud Functions entrypoint (minimal): reuse the main Express app
const app = require("./server");
exports.main = app;
