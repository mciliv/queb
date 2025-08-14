#!/usr/bin/env node

// Dev Browser Manager
// - Launches a dedicated Puppeteer browser for the app
// - Closes the browser as soon as source changes are detected
// - Reopens the browser when the frontend build output updates

const fs = require('fs');
const path = require('path');
const http = require('http');
const puppeteer = require('puppeteer');

const PROJECT_ROOT = path.resolve(__dirname, '..');
const FRONTEND_SRC_DIR = path.join(PROJECT_ROOT, 'frontend');
const FRONTEND_DIST_FILE = path.join(PROJECT_ROOT, 'frontend', 'dist', 'bundle.js');
const TARGET_URL = process.env.FRONTEND_URL || `http://localhost:${process.env.PORT || 3000}`;
const USER_DATA_DIR = path.join(PROJECT_ROOT, 'test', `chrome-molecular-profile-${Date.now()}-dev`);
const PID_FILE = '/tmp/dev_browser_pid';

let browser = null;
let page = null;
let closing = false;
let debounceTimer = null;

function log(msg) {
  console.log(msg);
}

function writePidFile() {
  try {
    fs.writeFileSync(PID_FILE, String(process.pid));
  } catch (_) {}
}

async function waitForServer(url, timeoutMs = 20000) {
  const deadline = Date.now() + timeoutMs;
  const urlObj = new URL(url);
  return new Promise((resolve, reject) => {
    const tryOnce = () => {
      const req = http.request({
        method: 'GET',
        hostname: urlObj.hostname,
        port: urlObj.port || 80,
        path: urlObj.pathname,
      }, (res) => {
        if (res.statusCode && res.statusCode >= 200 && res.statusCode < 500) {
          resolve(true);
        } else {
          next();
        }
      });
      req.on('error', next);
      req.end();
    };
    const next = () => {
      if (Date.now() > deadline) return reject(new Error('Server not ready'));
      setTimeout(tryOnce, 500);
    };
    tryOnce();
  });
}

async function openBrowser() {
  if (browser) return; // Already open
  await waitForServer(TARGET_URL).catch(() => {});
  browser = await puppeteer.launch({
    headless: false,
    defaultViewport: { width: 1600, height: 1000 },
    userDataDir: USER_DATA_DIR,
    args: [
      '--no-sandbox',
      '--disable-setuid-sandbox',
      '--disable-web-security',
      '--no-first-run',
      '--disable-default-apps',
      '--disable-infobars',
      '--disable-background-timer-throttling',
      '--disable-backgrounding-occluded-windows',
      '--disable-renderer-backgrounding'
    ]
  });
  page = await browser.newPage();
  await page.goto(TARGET_URL, { waitUntil: 'networkidle0' }).catch(() => {});
  log('🌐 Dev browser opened');
}

async function closeBrowser() {
  if (!browser || closing) return;
  closing = true;
  try {
    await browser.close();
  } catch (_) {}
  browser = null;
  page = null;
  closing = false;
  log('🧹 Dev browser closed');
}

function debounce(fn, delay) {
  return function debounced(...args) {
    if (debounceTimer) clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => fn.apply(this, args), delay);
  };
}

function watchSources() {
  try {
    fs.watch(FRONTEND_SRC_DIR, { recursive: true }, debounce(async () => {
      await closeBrowser();
    }, 150));
    log('👀 Watching frontend sources for changes');
  } catch (err) {
    log(`⚠️ Recursive watch not supported: ${err.message}`);
    // Fallback: watch key subdirectories individually
    const subdirs = ['core', 'components', 'assets'];
    subdirs.forEach((d) => {
      const dir = path.join(FRONTEND_SRC_DIR, d);
      if (fs.existsSync(dir)) {
        fs.watch(dir, { recursive: true }, debounce(async () => {
          await closeBrowser();
        }, 150));
      }
    });
  }
}

function watchBuildArtifact() {
  if (!fs.existsSync(FRONTEND_DIST_FILE)) {
    try { fs.mkdirSync(path.dirname(FRONTEND_DIST_FILE), { recursive: true }); } catch (_) {}
    fs.writeFileSync(FRONTEND_DIST_FILE, '');
  }
  fs.watchFile(FRONTEND_DIST_FILE, { interval: 200 }, debounce(async () => {
    // Reopen only if currently closed
    if (!browser) {
      await openBrowser();
    }
  }, 100));
  log('🔄 Watching build output for completion');
}

async function main() {
  writePidFile();
  watchSources();
  watchBuildArtifact();
  // Initial open once server is ready
  await openBrowser().catch(() => {});

  const cleanup = async () => {
    fs.unwatchFile(FRONTEND_DIST_FILE);
    await closeBrowser();
    try { fs.unlinkSync(PID_FILE); } catch (_) {}
    process.exit(0);
  };

  process.on('SIGINT', cleanup);
  process.on('SIGTERM', cleanup);
}

main().catch((err) => {
  log(`❌ Dev browser manager error: ${err.message}`);
  process.exit(1);
});





