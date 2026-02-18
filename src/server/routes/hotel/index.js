const express = require('express');
const fs = require('fs');
const path = require('path');

/**
 * Setup hotel card routes for the queb application
 * Integrates secure video hosting and hotel request card functionality
 */
function setupHotelRoutes(app, { logger }) {
  const HOTEL_LOG_FILE = path.join(__dirname, '..', '..', '..', '..', 'logs', 'hotel', 'access.log');
  const HOTEL_ADMIN_KEY_FILE = path.join(__dirname, '..', '..', '..', '..', '.hotel_admin_key');

  // Get admin key from file or environment
  let adminKey = process.env.HOTEL_ADMIN_KEY;
  if (!adminKey && fs.existsSync(HOTEL_ADMIN_KEY_FILE)) {
    try {
      adminKey = fs.readFileSync(HOTEL_ADMIN_KEY_FILE, 'utf8').trim();
    } catch (error) {
      logger.warn('Failed to read hotel admin key file:', error.message);
      adminKey = 'change-me-please'; // fallback
    }
  } else if (!adminKey) {
    adminKey = 'change-me-please'; // default fallback
  }

  // Security headers middleware for hotel routes
  const hotelSecurityHeaders = (req, res, next) => {
    // Prevent clickjacking
    res.setHeader('X-Frame-Options', 'DENY');

    // Prevent MIME sniffing
    res.setHeader('X-Content-Type-Options', 'nosniff');

    // XSS protection
    res.setHeader('X-XSS-Protection', '1; mode=block');

    // Content Security Policy - only allow same origin
    res.setHeader('Content-Security-Policy', "default-src 'self'; media-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'");

    // Disable caching for sensitive content
    res.setHeader('Cache-Control', 'no-store, no-cache, must-revalidate, private');
    res.setHeader('Pragma', 'no-cache');
    res.setHeader('Expires', '0');

    // Remove server fingerprinting
    res.removeHeader('X-Powered-By');

    next();
  };

  // Apply security headers to all hotel routes
  app.use('/hotel', hotelSecurityHeaders);

  // Viewer tracking middleware
  function trackHotelAccess(resourceType) {
    return (req, res, next) => {
      const timestamp = new Date().toISOString();
      const ip = req.headers['x-forwarded-for'] || req.socket.remoteAddress;
      const userAgent = req.headers['user-agent'] || 'Unknown';
      const referer = req.headers['referer'] || 'Direct';

      const logEntry = {
        timestamp,
        resource: resourceType,
        ip,
        userAgent,
        referer,
        path: req.path,
        module: 'hotel'
      };

      const logLine = JSON.stringify(logEntry) + '\n';

      fs.appendFile(HOTEL_LOG_FILE, logLine, (err) => {
        if (err) logger.error('Failed to write hotel access log:', err);
      });

      next();
    };
  }

  // Ensure hotel logs directory exists
  const hotelLogsDir = path.dirname(HOTEL_LOG_FILE);
  if (!fs.existsSync(hotelLogsDir)) {
    fs.mkdirSync(hotelLogsDir, { recursive: true });
  }

  // Hotel card routes
  app.get('/hotel', trackHotelAccess('hotel-card'), (req, res) => {
    const hotelCardPath = path.join(__dirname, '..', '..', '..', '..', 'public', 'hotel', 'hotel-card.html');
    res.sendFile(hotelCardPath, (err) => {
      if (err) {
        logger.error('Failed to serve hotel card:', err);
        res.status(500).send('Hotel card not available');
      }
    });
  });

  // Hotel logs viewer (secured endpoint)
  app.get('/hotel/admin/logs', (req, res) => {
    const requestKey = req.query.key;

    if (requestKey !== adminKey) {
      return res.status(403).send('Access denied. Use ?key=YOUR_HOTEL_ADMIN_KEY');
    }

    const logsHtmlPath = path.join(__dirname, '..', '..', '..', '..', 'public', 'hotel', 'logs.html');
    res.sendFile(logsHtmlPath, (err) => {
      if (err) {
        logger.error('Failed to serve hotel logs:', err);
        res.status(500).send('Logs dashboard not available');
      }
    });
  });

  // API endpoint to fetch hotel logs data
  app.get('/hotel/admin/api/logs', (req, res) => {
    const requestKey = req.query.key;

    if (requestKey !== adminKey) {
      return res.status(403).json({ error: 'Access denied' });
    }

    if (!fs.existsSync(HOTEL_LOG_FILE)) {
      return res.json([]);
    }

    fs.readFile(HOTEL_LOG_FILE, 'utf8', (err, data) => {
      if (err) {
        logger.error('Failed to read hotel logs:', err);
        return res.status(500).json({ error: 'Failed to read logs' });
      }

      const logs = data
        .trim()
        .split('\n')
        .filter(line => line)
        .map(line => {
          try {
            return JSON.parse(line);
          } catch (e) {
            return null;
          }
        })
        .filter(log => log !== null)
        .reverse(); // Most recent first

      res.json(logs);
    });
  });

  // Video serving with security (if video exists)
  app.get('/hotel/video/:filename', trackHotelAccess('hotel-video'), (req, res) => {
    const filename = req.params.filename;
    const videoPath = path.join(__dirname, '..', '..', '..', '..', 'public', 'hotel', filename);

    if (!fs.existsSync(videoPath)) {
      return res.status(404).send('Video not found');
    }

    // Additional security for video files
    res.setHeader('Content-Type', 'video/mp4');
    res.setHeader('Accept-Ranges', 'bytes');

    res.sendFile(videoPath, (err) => {
      if (err) {
        logger.error('Failed to serve hotel video:', err);
        res.status(500).send('Video unavailable');
      }
    });
  });

  logger.info('Hotel routes integrated successfully');
}

module.exports = { setupHotelRoutes };