/**
 * User management routes
 * Handles user creation, updates, and retrieval by device token
 */

/**
 * Setup user management routes
 * @param {express.Application} app - Express application instance
 * @param {ServiceContainer} container - Dependency injection container
 */
async function setupUserRoutes(app, container) {
  const userService = await container.get('userService');
  const logger = await container.get('logger');
  
  // Create/update user
  app.post('/api/user', async (req, res, next) => {
    try {
      const { deviceToken, ...userData } = req.body;
      
      if (!deviceToken) {
        return res.status(400).json({ 
          error: 'Device token required' 
        });
      }
      
      let user = await userService.getUserByDeviceToken(deviceToken);
      
      if (user) {
        // Update existing user
        await userService.updateUser(deviceToken, userData);
        user = await userService.getUserByDeviceToken(deviceToken);
        logger.info('User updated', { userId: user.id });
      } else {
        user = await userService.createUser({ deviceToken, ...userData });
        logger.info('User created', { userId: user.id });
      }
      
      res.json(user);
    } catch (error) {
      next(error);
    }
  });
  
  // Get user
  app.get('/api/user/:deviceToken', async (req, res, next) => {
    try {
      const { deviceToken } = req.params;
      const user = await userService.getUserByDeviceToken(deviceToken);
      
      if (!user) {
        return res.status(404).json({ 
          error: 'User not found' 
        });
      }
      
      res.json(user);
    } catch (error) {
      next(error);
    }
  });
}

module.exports = {
  setupUserRoutes
};






