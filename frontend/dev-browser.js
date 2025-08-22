#!/usr/bin/env node

// Dev Browser Manager
// - Launches a dedicated Puppeteer browser for the app
// - Closes the browser as soon as source changes are detected
// - Reopens the browser when the frontend build output updates
// - In human verification mode, do not auto-close; instead reload the page

const fs = require('fs');
const path = require('path');
const http = require('http');
const { exec } = require('child_process');
const puppeteer = require('puppeteer');

const PROJECT_ROOT = path.resolve(__dirname, '..');
const FRONTEND_SRC_DIR = path.join(PROJECT_ROOT, 'frontend');
const FRONTEND_DIST_FILE = path.join(PROJECT_ROOT, 'frontend', 'dist', 'bundle.js');
const TARGET_URL = process.env.FRONTEND_URL || `http://localhost:${process.env.PORT || 3000}`;
const USER_DATA_DIR = path.join(PROJECT_ROOT, 'test', `chrome-molecular-profile-${Date.now()}-dev`);
const PID_FILE = '/tmp/dev_browser_pid';
const STICKY = (
  process.env.VISUAL_STICKY === '1' ||
  process.env.HUMAN_VERIFY === '1' ||
  process.env.VISUAL_MODE === 'on' ||
  process.env.DEV_VISUAL === 'on'
);

let browser = null;
let page = null;
let closing = false;
let debounceTimer = null;
let usedSystemBrowserFallback = false;

function removeUserDataDir() {
  try {
    if (fs.existsSync(USER_DATA_DIR)) {
      fs.rmSync(USER_DATA_DIR, { recursive: true, force: true });
    }
  } catch (_) {}
}

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
  try {
    await waitForServer(TARGET_URL, 60000);
  } catch (_) {
    log('⏳ Server not ready, postponing browser launch');
    return;
  }
  try {
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
  } catch (err) {
    if (!usedSystemBrowserFallback) {
      const opener = process.platform === 'darwin' ? 'open' : (process.platform === 'win32' ? 'start' : 'xdg-open');
      try {
        exec(`${opener} "${TARGET_URL}"`, { stdio: 'ignore' });
        usedSystemBrowserFallback = true;
        log('🌐 Dev browser opened (system)');
      } catch (_) {
        log(`⚠️ Failed to launch Puppeteer and system browser`);
      }
    }
    return;
  }
  const pages = await browser.pages();
  page = pages && pages.length > 0 ? pages[0] : await browser.newPage();
  // Close any extra initial tabs (often about:blank)
  if (pages && pages.length > 1) {
    for (const extra of pages.slice(1)) {
      try { await extra.close(); } catch (_) {}
    }
  }
  await page.bringToFront().catch(() => {});
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
    // Ignore source changes to avoid double-reload; rely on build artifact watcher
    fs.watch(FRONTEND_SRC_DIR, { recursive: true }, debounce(async () => {
      // No action; build watcher will handle reload
    }, 500));
    log('👀 Ignoring source changes (build watcher will trigger reload)');
  } catch (err) {
    log(`⚠️ Recursive watch not supported: ${err.message}`);
    // Fallback: watch key subdirectories individually
    const subdirs = ['core', 'components', 'assets'];
    subdirs.forEach((d) => {
      const dir = path.join(FRONTEND_SRC_DIR, d);
      if (fs.existsSync(dir)) {
        fs.watch(dir, { recursive: true }, debounce(async () => {
          // No action; build watcher will handle reload
        }, 500));
      }
    });
  }
}

function watchBuildArtifact() {
  if (!fs.existsSync(FRONTEND_DIST_FILE)) {
    try { fs.mkdirSync(path.dirname(FRONTEND_DIST_FILE), { recursive: true }); } catch (_) {}
    fs.writeFileSync(FRONTEND_DIST_FILE, '');
  }
  fs.watchFile(FRONTEND_DIST_FILE, { interval: 250 }, debounce(async () => {
    if (STICKY && browser && page) {
      try {
        await page.reload({ waitUntil: 'networkidle0' });
        await page.bringToFront().catch(() => {});
        log('🔁 Dev browser reloaded after build (sticky)');
      } catch (_) {}
    } else {
      // Reopen only if currently closed
      if (!browser) {
        await openBrowser();
      }
    }
  }, 400));
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
    // Remove the temporary Chrome profile used for this dev run
    removeUserDataDir();
    process.exit(0);
  };

  process.on('SIGINT', cleanup);
  process.on('SIGTERM', cleanup);
  // Best-effort sync cleanup on process exit
  process.on('exit', () => {
    try { fs.unwatchFile(FRONTEND_DIST_FILE); } catch (_) {}
    try { fs.unlinkSync(PID_FILE); } catch (_) {}
    removeUserDataDir();
  });
}

main().catch((err) => {
  log(`❌ Dev browser manager error: ${err.message}`);
  process.exit(1);
});





