// Backend Node.js Utilities
// =========================

const fs = require('fs').promises;
const path = require('path');

// Project Setup Utilities
const ProjectSetup = {
    // Initialize Node.js project
    initProject: async (projectName, options = {}) => {
        const packageJson = {
            name: projectName,
            version: "1.0.0",
            description: options.description || "",
            main: options.main || "index.js",
            scripts: {
                start: "node index.js",
                dev: "nodemon index.js",
                test: "jest",
                ...options.scripts
            },
            dependencies: {
                express: "^4.18.2",
                cors: "^2.8.5",
                ...options.dependencies
            },
            devDependencies: {
                nodemon: "^3.0.0",
                jest: "^29.0.0",
                ...options.devDependencies
            }
        };

        await fs.writeFile('package.json', JSON.stringify(packageJson, null, 2));
        console.log('✅ package.json created');
    },

    // Install dependencies
    installDeps: async () => {
        const { exec } = require('child_process');
        return new Promise((resolve, reject) => {
            exec('npm install', (error, stdout, stderr) => {
                if (error) {
                    console.error('❌ Failed to install dependencies:', error);
                    reject(error);
                } else {
                    console.log('✅ Dependencies installed');
                    resolve(stdout);
                }
            });
        });
    },

    // Run tests
    runTests: async () => {
        const { exec } = require('child_process');
        return new Promise((resolve, reject) => {
            exec('npm test', (error, stdout, stderr) => {
                if (error) {
                    console.error('❌ Tests failed:', error);
                    reject(error);
                } else {
                    console.log('✅ Tests passed');
                    resolve(stdout);
                }
            });
        });
    }
};

// Express Utilities
const ExpressUtils = {
    // Create basic Express app
    createApp: (options = {}) => {
        const express = require('express');
        const cors = require('cors');
        
        const app = express();
        
        // Middleware
        app.use(cors());
        app.use(express.json());
        app.use(express.urlencoded({ extended: true }));
        
        // Custom middleware
        if (options.middleware) {
            options.middleware.forEach(mw => app.use(mw));
        }
        
        // Error handling
        app.use((err, req, res, next) => {
            console.error(err.stack);
            res.status(500).json({ error: 'Something went wrong!' });
        });
        
        return app;
    },

    // Add routes from directory
    addRoutesFromDir: async (app, routesDir) => {
        try {
            const files = await fs.readdir(routesDir);
            for (const file of files) {
                if (file.endsWith('.js')) {
                    const route = require(path.join(routesDir, file));
                    app.use(route);
                }
            }
        } catch (error) {
            console.error('Failed to load routes:', error);
        }
    },

    // Health check endpoint
    addHealthCheck: (app) => {
        app.get('/health', (req, res) => {
            res.json({ 
                status: 'ok', 
                timestamp: new Date().toISOString(),
                uptime: process.uptime()
            });
        });
    }
};

// Database Utilities
const DatabaseUtils = {
    // Simple in-memory database
    createMemoryDB: () => {
        const db = new Map();
        
        return {
            set: (key, value) => db.set(key, value),
            get: (key) => db.get(key),
            delete: (key) => db.delete(key),
            has: (key) => db.has(key),
            clear: () => db.clear(),
            size: () => db.size,
            keys: () => Array.from(db.keys()),
            values: () => Array.from(db.values()),
            entries: () => Array.from(db.entries())
        };
    },

    // File-based database
    createFileDB: async (filePath) => {
        let data = {};
        
        try {
            const content = await fs.readFile(filePath, 'utf8');
            data = JSON.parse(content);
        } catch (error) {
            // File doesn't exist, start with empty data
        }
        
        return {
            set: async (key, value) => {
                data[key] = value;
                await fs.writeFile(filePath, JSON.stringify(data, null, 2));
            },
            get: (key) => data[key],
            delete: async (key) => {
                delete data[key];
                await fs.writeFile(filePath, JSON.stringify(data, null, 2));
            },
            has: (key) => key in data,
            clear: async () => {
                data = {};
                await fs.writeFile(filePath, JSON.stringify(data, null, 2));
            },
            keys: () => Object.keys(data),
            values: () => Object.values(data),
            entries: () => Object.entries(data)
        };
    }
};

// Validation Utilities
const Validation = {
    // Email validation
    isEmail: (email) => {
        const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
        return emailRegex.test(email);
    },

    // URL validation
    isURL: (url) => {
        try {
            new URL(url);
            return true;
        } catch {
            return false;
        }
    },

    // Required field validation
    isRequired: (value) => {
        return value !== null && value !== undefined && value.toString().trim() !== '';
    },

    // Validate request body
    validateBody: (schema) => {
        return (req, res, next) => {
            const errors = [];
            
            for (const [field, rules] of Object.entries(schema)) {
                const value = req.body[field];
                
                if (rules.required && !Validation.isRequired(value)) {
                    errors.push(`${field} is required`);
                }
                
                if (value && rules.type === 'email' && !Validation.isEmail(value)) {
                    errors.push(`${field} must be a valid email`);
                }
                
                if (value && rules.type === 'url' && !Validation.isURL(value)) {
                    errors.push(`${field} must be a valid URL`);
                }
            }
            
            if (errors.length > 0) {
                return res.status(400).json({ errors });
            }
            
            next();
        };
    }
};

// Logging Utilities
const Logger = {
    info: (message, data = {}) => {
        console.log(`[INFO] ${new Date().toISOString()} - ${message}`, data);
    },
    
    error: (message, error = {}) => {
        console.error(`[ERROR] ${new Date().toISOString()} - ${message}`, error);
    },
    
    warn: (message, data = {}) => {
        console.warn(`[WARN] ${new Date().toISOString()} - ${message}`, data);
    },
    
    debug: (message, data = {}) => {
        if (process.env.NODE_ENV === 'development') {
            console.log(`[DEBUG] ${new Date().toISOString()} - ${message}`, data);
        }
    }
};

// Security Utilities
const Security = {
    // Generate random string
    randomString: (length = 32) => {
        const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
        let result = '';
        for (let i = 0; i < length; i++) {
            result += chars.charAt(Math.floor(Math.random() * chars.length));
        }
        return result;
    },

    // Basic rate limiting
    rateLimit: (maxRequests = 100, windowMs = 15 * 60 * 1000) => {
        const requests = new Map();
        
        return (req, res, next) => {
            const ip = req.ip;
            const now = Date.now();
            const windowStart = now - windowMs;
            
            // Clean old requests
            if (requests.has(ip)) {
                requests.set(ip, requests.get(ip).filter(time => time > windowStart));
            } else {
                requests.set(ip, []);
            }
            
            const currentRequests = requests.get(ip);
            
            if (currentRequests.length >= maxRequests) {
                return res.status(429).json({ error: 'Too many requests' });
            }
            
            currentRequests.push(now);
            next();
        };
    }
};

module.exports = {
    ProjectSetup,
    ExpressUtils,
    DatabaseUtils,
    Validation,
    Logger,
    Security
}; 