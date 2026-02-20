const express = require('express');
const fs = require('fs');
const path = require('path');

const app = express();
const PORT = process.env.PORT || 3333;
const LOG_FILE = path.join(__dirname, 'logs', 'access.log');

// Security headers
app.use((req, res, next) => {
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
});

// Viewer tracking middleware
function trackAccess(resourceType) {
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
      path: req.path
    };

    const logLine = JSON.stringify(logEntry) + '\n';

    fs.appendFile(LOG_FILE, logLine, (err) => {
      if (err) console.error('Failed to write log:', err);
    });

    next();
  };
}

// Serve static files from public directory
app.use(express.static('public'));

// Video access tracking
app.get('/video/cc472.mp4', trackAccess('video'), (req, res) => {
  const videoPath = path.join(__dirname, 'public', 'cc472.mp4');

  if (!fs.existsSync(videoPath)) {
    return res.status(404).send('Video not found');
  }

  res.sendFile(videoPath);
});

// Hotel Request Card - main page
app.get('/', trackAccess('hotel-card'), (req, res) => {
  res.sendFile(path.join(__dirname, 'public', 'hotel-card.html'));
});

// Also serve at /hotel for clarity
app.get('/hotel', trackAccess('hotel-card'), (req, res) => {
  res.sendFile(path.join(__dirname, 'public', 'hotel-card.html'));
});

// Logs viewer (secured endpoint)
app.get('/admin/logs', (req, res) => {
  // Basic auth check - you should set ADMIN_KEY environment variable
  const adminKey = req.query.key;
  const expectedKey = process.env.ADMIN_KEY || 'change-me-please';

  if (adminKey !== expectedKey) {
    return res.status(403).send('Access denied. Use ?key=YOUR_ADMIN_KEY');
  }

  res.sendFile(path.join(__dirname, 'public', 'logs.html'));
});

// API endpoint to fetch logs data
app.get('/admin/api/logs', (req, res) => {
  const adminKey = req.query.key;
  const expectedKey = process.env.ADMIN_KEY || 'change-me-please';

  if (adminKey !== expectedKey) {
    return res.status(403).json({ error: 'Access denied' });
  }

  if (!fs.existsSync(LOG_FILE)) {
    return res.json([]);
  }

  fs.readFile(LOG_FILE, 'utf8', (err, data) => {
    if (err) {
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

// Ensure logs directory exists
if (!fs.existsSync(path.join(__dirname, 'logs'))) {
  fs.mkdirSync(path.join(__dirname, 'logs'));
}

app.listen(PORT, '0.0.0.0', () => {
  console.log(`\n🏨 Hotel Request Card running at:`);
  console.log(`   Local:    http://localhost:${PORT}`);
  console.log(`   Network:  http://0.0.0.0:${PORT}`);
  console.log(`   Domain:   https://queb.space (if tunneled)`);
  console.log(`\n📱 Phone can access at: https://queb.space/hotel`);
  console.log(`\n✓ No-cache headers enabled - phone always gets latest version`);
});
