/**
 * MCP Integration Service
 * 
 * Provides utility functions to leverage MCP-like capabilities in the application.
 * While MCPs are primarily for AI assistants, this service provides similar
 * functionality for runtime operations.
 * 
 * Key capabilities:
 * - Documentation lookup for chemical libraries
 * - Enhanced database operations with schema introspection
 * - Structured file operations
 */

const { Pool } = require('pg');

class MCPIntegrationService {
  constructor(config = {}) {
    this.dbPool = config.dbPool || null;
    this.logger = config.logger || console;
  }

  /**
   * Initialize database connection pool if not provided
   */
  initializeDatabase() {
    if (!this.dbPool && process.env.DATABASE_URL) {
      this.dbPool = new Pool({
        connectionString: process.env.DATABASE_URL,
        max: 10,
        idleTimeoutMillis: 30000,
      });
      this.logger.info('Database pool initialized for MCP service');
    }
  }

  /**
   * Get database schema information
   * Useful for understanding data structure and validation
   */
  async getDatabaseSchema(tableName = null) {
    this.initializeDatabase();
    
    if (!this.dbPool) {
      throw new Error('Database not configured');
    }

    try {
      const query = tableName
        ? `
          SELECT 
            column_name, 
            data_type, 
            is_nullable,
            column_default
          FROM information_schema.columns
          WHERE table_name = $1
          ORDER BY ordinal_position;
        `
        : `
          SELECT 
            table_name,
            table_type
          FROM information_schema.tables
          WHERE table_schema = 'public'
          ORDER BY table_name;
        `;

      const params = tableName ? [tableName] : [];
      const result = await this.dbPool.query(query, params);
      
      return result.rows;
    } catch (error) {
      this.logger.error('Failed to get database schema:', error);
      throw error;
    }
  }

  /**
   * Execute safe read-only database queries
   * Prevents accidental modifications
   */
  async executeReadQuery(sql, params = []) {
    this.initializeDatabase();
    
    if (!this.dbPool) {
      throw new Error('Database not configured');
    }

    // Basic safety check - reject queries that look like modifications
    const upperSQL = sql.trim().toUpperCase();
    const dangerousKeywords = ['INSERT', 'UPDATE', 'DELETE', 'DROP', 'ALTER', 'CREATE', 'TRUNCATE'];
    
    if (dangerousKeywords.some(keyword => upperSQL.startsWith(keyword))) {
      throw new Error('Only SELECT queries are allowed in executeReadQuery');
    }

    try {
      const result = await this.dbPool.query(sql, params);
      return result.rows;
    } catch (error) {
      this.logger.error('Query execution failed:', error);
      throw error;
    }
  }

  /**
   * Get table statistics
   * Useful for monitoring and optimization
   */
  async getTableStats(tableName) {
    this.initializeDatabase();
    
    if (!this.dbPool) {
      throw new Error('Database not configured');
    }

    try {
      const query = `
        SELECT 
          schemaname,
          relname as table_name,
          n_live_tup as row_count,
          n_dead_tup as dead_rows,
          last_vacuum,
          last_autovacuum,
          last_analyze,
          last_autoanalyze
        FROM pg_stat_user_tables
        WHERE relname = $1;
      `;
      
      const result = await this.dbPool.query(query, [tableName]);
      return result.rows[0] || null;
    } catch (error) {
      this.logger.error('Failed to get table stats:', error);
      throw error;
    }
  }

  /**
   * Search for molecular data patterns in database
   * Enhanced search with chemical-specific logic
   */
  async searchMolecularData(searchTerm, options = {}) {
    this.initializeDatabase();
    
    if (!this.dbPool) {
      throw new Error('Database not configured');
    }

    const limit = options.limit || 50;
    const offset = options.offset || 0;

    try {
      // This is a template - adjust based on actual schema
      const query = `
        SELECT *
        FROM molecules
        WHERE 
          name ILIKE $1 OR
          smiles ILIKE $1 OR
          formula ILIKE $1
        ORDER BY created_at DESC
        LIMIT $2 OFFSET $3;
      `;
      
      const result = await this.dbPool.query(query, [`%${searchTerm}%`, limit, offset]);
      return result.rows;
    } catch (error) {
      // Table might not exist yet
      if (error.code === '42P01') {
        this.logger.warn('Molecules table does not exist yet');
        return [];
      }
      throw error;
    }
  }

  /**
   * Enhanced file operations with validation
   * Provides safer file handling than direct fs operations
   */
  async safeReadFile(filePath, options = {}) {
    const fs = require('fs').promises;
    const path = require('path');
    
    // Validate file path is within project
    const projectRoot = process.cwd();
    const absolutePath = path.resolve(projectRoot, filePath);
    
    if (!absolutePath.startsWith(projectRoot)) {
      throw new Error('File path must be within project directory');
    }

    try {
      const content = await fs.readFile(absolutePath, options.encoding || 'utf8');
      return {
        path: filePath,
        absolutePath,
        content,
        size: Buffer.byteLength(content),
      };
    } catch (error) {
      this.logger.error('Failed to read file:', error);
      throw error;
    }
  }

  /**
   * Enhanced file writing with backup
   */
  async safeWriteFile(filePath, content, options = {}) {
    const fs = require('fs').promises;
    const path = require('path');
    
    // Validate file path is within project
    const projectRoot = process.cwd();
    const absolutePath = path.resolve(projectRoot, filePath);
    
    if (!absolutePath.startsWith(projectRoot)) {
      throw new Error('File path must be within project directory');
    }

    try {
      // Create backup if file exists
      if (options.backup) {
        try {
          await fs.access(absolutePath);
          const backupPath = `${absolutePath}.backup`;
          await fs.copyFile(absolutePath, backupPath);
          this.logger.info(`Created backup at ${backupPath}`);
        } catch (err) {
          // File doesn't exist, no backup needed
        }
      }

      await fs.writeFile(absolutePath, content, options.encoding || 'utf8');
      
      return {
        path: filePath,
        absolutePath,
        size: Buffer.byteLength(content),
        success: true,
      };
    } catch (error) {
      this.logger.error('Failed to write file:', error);
      throw error;
    }
  }

  /**
   * List SDF files in a directory
   * Chemical structure file management
   */
  async listSDFFiles(directory = 'src/backend/sdf_files') {
    const fs = require('fs').promises;
    const path = require('path');
    
    const projectRoot = process.cwd();
    const absolutePath = path.resolve(projectRoot, directory);
    
    if (!absolutePath.startsWith(projectRoot)) {
      throw new Error('Directory must be within project');
    }

    try {
      const files = await fs.readdir(absolutePath);
      const sdfFiles = files.filter(file => file.endsWith('.sdf'));
      
      const fileDetails = await Promise.all(
        sdfFiles.map(async (file) => {
          const filePath = path.join(absolutePath, file);
          const stats = await fs.stat(filePath);
          
          return {
            name: file,
            path: path.join(directory, file),
            size: stats.size,
            modified: stats.mtime,
            moleculeName: file.replace('.sdf', ''),
          };
        })
      );
      
      return fileDetails;
    } catch (error) {
      if (error.code === 'ENOENT') {
        this.logger.warn(`Directory ${directory} does not exist`);
        return [];
      }
      throw error;
    }
  }

  /**
   * Cleanup method
   */
  async cleanup() {
    if (this.dbPool) {
      await this.dbPool.end();
      this.logger.info('Database pool closed');
    }
  }
}

// Singleton instance
let instance = null;

function getMCPService(config = {}) {
  if (!instance) {
    instance = new MCPIntegrationService(config);
  }
  return instance;
}

module.exports = {
  MCPIntegrationService,
  getMCPService,
};









