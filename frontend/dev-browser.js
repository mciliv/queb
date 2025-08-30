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
const TARGET_URL = process.env.FRONTEND_URL || `https://localhost:${process.env.PORT || 3001}`;
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
      // Open DevTools automatically
      devtools: true,
      // Use full screen window size
      defaultViewport: null,
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
        '--disable-renderer-backgrounding',
        '--start-maximized',
        '--window-size=1920,1080',
        '--ignore-certificate-errors',
        '--ignore-ssl-errors',
        // Auto-open DevTools for tabs
        '--auto-open-devtools-for-tabs'
      ]
    });
  } catch (err) {
    if (!usedSystemBrowserFallback) {
      try {
        if (process.platform === 'darwin') {
          // Try to open Google Chrome in maximized window with DevTools
          exec(
            `open -a "Google Chrome" --args --start-maximized --window-size=1920,1080 --ignore-certificate-errors --ignore-ssl-errors --auto-open-devtools-for-tabs "${TARGET_URL}"`,
            { stdio: 'ignore' }
          );
        } else if (process.platform === 'win32') {
          exec(
            `start chrome --start-maximized --window-size=1920,1080 --ignore-certificate-errors --ignore-ssl-errors --auto-open-devtools-for-tabs "${TARGET_URL}"`,
            { shell: true, stdio: 'ignore' }
          );
        } else {
          // Fallback: open default browser
          exec(`xdg-open "${TARGET_URL}"`, { stdio: 'ignore' });
        }
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
  // Set window to fill screen without going fullscreen
  try {
    await page.setViewport({ width: 1920, height: 1080 });
  } catch (_) {}
  // Ensure body/root consume full viewport
  try {
    await page.addStyleTag({ content: 'html,body,#root{height:100vh;width:100vw;margin:0;padding:0;overflow:hidden;}' });
  } catch (_) {}
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
    if (browser && page) {
      try {
        // Always reload the tab when build artifact changes
        await page.reload({ waitUntil: 'networkidle0' });
        await page.bringToFront().catch(() => {});
        log('🔁 Dev browser reloaded after build');
      } catch (_) {}
    } else {
      // Reopen only if currently closed
      await openBrowser();
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
