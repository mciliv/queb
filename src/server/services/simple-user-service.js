// simple-user-service.js - In-memory user management for demo
class SimpleUserService {
  constructor() {
    this.users = new Map();
    this.usageStats = new Map();

  }

  async initializeTables() {
    // No-op for in-memory storage

    return true;
  }

  async createUser(userData) {
    const { deviceToken, deviceInfo, name } = userData;
    
    if (this.users.has(deviceToken)) {
      throw new Error('Device token already exists');
    }

    const user = {
      id: this.users.size + 1,
      device_token: deviceToken,
      device_info: deviceInfo,
      name: name || null,
      usage: 0,
      created_at: new Date().toISOString(),
      last_used: new Date().toISOString(),
      is_active: true
    };

    this.users.set(deviceToken, user);

    return user;
  }

  async getUserByDeviceToken(deviceToken) {
    const user = this.users.get(deviceToken);
    if (user) {
      // Update last_used
      user.last_used = new Date().toISOString();
    }
    return user || null;
  }

  async updateUser(deviceToken, updateData) {
    const user = this.users.get(deviceToken);
    if (!user) {
      throw new Error('User not found');
    }

    // Update fields
    Object.assign(user, updateData);
    user.last_used = new Date().toISOString();
    
    this.users.set(deviceToken, user);
    return user;
  }

  async incrementUsage(deviceToken, amount = 1) {
    const user = this.users.get(deviceToken);
    if (user) {
      user.usage += amount;
      user.last_used = new Date().toISOString();
  
    }
    return user;
  }

  async getAllUsers() {
    return Array.from(this.users.values());
  }

  async deleteUser(deviceToken) {
    const deleted = this.users.delete(deviceToken);
    if (deleted) {
  
    }
    return deleted;
  }

  // Stats for monitoring
  getStats() {
    const users = Array.from(this.users.values());
    return {
      totalUsers: users.length,
      activeUsers: users.filter(u => u.is_active).length,
      totalUsage: users.reduce((sum, u) => sum + u.usage, 0),
      averageUsage: users.length > 0 ? users.reduce((sum, u) => sum + u.usage, 0) / users.length : 0
    };
  }
}

module.exports = SimpleUserService; 