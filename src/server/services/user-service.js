// user-service.js - PostgreSQL-based user management service

class UserService {
  constructor(pool) {
    this.pool = pool;
  }

  // Initialize database tables
  async initializeTables() {
    const client = await this.pool.connect();
    try {
      // Create users table if it doesn't exist
      await client.query(`
        CREATE TABLE IF NOT EXISTS users (
          id SERIAL PRIMARY KEY,
          device_token VARCHAR(255) UNIQUE NOT NULL,
          payment_method_id VARCHAR(255),
          device_info JSONB,
          name VARCHAR(255),
          usage INTEGER DEFAULT 0,
          created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
          last_used TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
          is_active BOOLEAN DEFAULT true
        );
      `);

      // Create indexes for performance
      await client.query(`
        CREATE INDEX IF NOT EXISTS idx_users_device_token ON users (device_token);
        CREATE INDEX IF NOT EXISTS idx_users_last_used ON users (last_used);
        CREATE INDEX IF NOT EXISTS idx_users_active ON users (is_active);
      `);


      return true;
    } catch (error) {

      throw error;
    } finally {
      client.release();
    }
  }

  // Create a new user
  async createUser(userData) {
    const client = await this.pool.connect();
    try {
      const {
        deviceToken,
        paymentMethodId,
        deviceInfo,
        name
      } = userData;

      const result = await client.query(`
        INSERT INTO users (device_token, payment_method_id, device_info, name)
        VALUES ($1, $2, $3, $4)
        RETURNING *;
      `, [deviceToken, paymentMethodId, deviceInfo, name]);


      return result.rows[0];
    } catch (error) {
      if (error.code === '23505') { // Unique constraint violation
        throw new Error('Device token already exists');
      }
      
      throw error;
    } finally {
      client.release();
    }
  }

  // Get user by device token
  async getUserByDeviceToken(deviceToken) {
    const client = await this.pool.connect();
    try {
      const result = await client.query(`
        SELECT * FROM users 
        WHERE device_token = $1 AND is_active = true;
      `, [deviceToken]);

      if (result.rows.length === 0) {
        return null;
      }

      // Update last_used timestamp
      await client.query(`
        UPDATE users 
        SET last_used = NOW() 
        WHERE device_token = $1;
      `, [deviceToken]);

      return result.rows[0];
    } catch (error) {
      
      throw error;
    } finally {
      client.release();
    }
  }

  // Update user information
  async updateUser(deviceToken, updateData) {
    const client = await this.pool.connect();
    try {
      const {
        paymentMethodId,
        deviceInfo,
        name
      } = updateData;

      const updates = [];
      const values = [];
      let paramCount = 0;

      if (paymentMethodId !== undefined) {
        updates.push(`payment_method_id = $${++paramCount}`);
        values.push(paymentMethodId);
      }
      
      if (deviceInfo !== undefined) {
        updates.push(`device_info = $${++paramCount}`);
        values.push(deviceInfo);
      }
      
      if (name !== undefined) {
        updates.push(`name = $${++paramCount}`);
        values.push(name);
      }

      if (updates.length === 0) {
        throw new Error('No update data provided');
      }

      // Always update last_used
      updates.push(`last_used = NOW()`);
      values.push(deviceToken);

      const result = await client.query(`
        UPDATE users 
        SET ${updates.join(', ')}
        WHERE device_token = $${++paramCount} AND is_active = true
        RETURNING *;
      `, values);

      if (result.rows.length === 0) {
        throw new Error('User not found');
      }

      const user = result.rows[0];

      return user;
    } catch (error) {

      throw error;
    } finally {
      client.release();
    }
  }

  // Increment usage counter for a user
  async incrementUsage(deviceToken) {
    const client = await this.pool.connect();
    try {
      const result = await client.query(`
        UPDATE users 
        SET usage = usage + 1, last_used = NOW()
        WHERE device_token = $1 AND is_active = true
        RETURNING usage, name;
      `, [deviceToken]);

      if (result.rows.length === 0) {
        throw new Error('User not found');
      }

      const user = result.rows[0];

      
      return user.usage;
    } catch (error) {

      throw error;
    } finally {
      client.release();
    }
  }

  // Get user statistics (for admin/debugging)
  async getUserStats() {
    const client = await this.pool.connect();
    try {
      const result = await client.query(`
        SELECT 
          COUNT(*) as total_users,
          SUM(usage) as total_analyses,
          AVG(usage) as avg_analyses_per_user,
          MAX(last_used) as last_activity
        FROM users 
        WHERE is_active = true;
      `);

      return result.rows[0];
    } catch (error) {

      throw error;
    } finally {
      client.release();
    }
  }

  // Cleanup old inactive users (for maintenance)
  async cleanupInactiveUsers(daysInactive = 90) {
    const client = await this.pool.connect();
    try {
      const result = await client.query(`
        UPDATE users 
        SET is_active = false
        WHERE last_used < NOW() - INTERVAL '${daysInactive} days'
        AND is_active = true
        RETURNING COUNT(*);
      `);

      const count = result.rowCount;
      if (count > 0) {
  
      }
      
      return count;
    } catch (error) {

      throw error;
    } finally {
      client.release();
    }
  }

  // Health check - verify database connection
  async healthCheck() {
    const client = await this.pool.connect();
    try {
      await client.query('SELECT 1');
      return { status: 'healthy', timestamp: new Date() };
    } catch (error) {
      return { status: 'unhealthy', error: error.message, timestamp: new Date() };
    } finally {
      client.release();
    }
  }
}

module.exports = UserService; 