// Stub: database recommender service
// Prompt: app.js line 439 requires this at startup via setupDatabaseRecommendationRoutes
class DatabaseRecommender {
  constructor({ logger, aiService, nameResolver }) {
    this.logger = logger;
    this.aiService = aiService;
    this.nameResolver = nameResolver;
  }

  async findChemicalContents(item) {
    throw new Error('DatabaseRecommender not yet implemented');
  }
}

module.exports = DatabaseRecommender;
