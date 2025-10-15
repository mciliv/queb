// Service Worker for Molecular Analysis PWA
const CACHE_NAME = 'mol-analysis-v1';
const STATIC_CACHE = 'mol-analysis-static-v1';

// Files to cache immediately
const STATIC_FILES = [
  '/',
  '/manifest.json',
  '/assets/style.css',
  '/assets/favicon.svg',
  '/dist/bundle.js',
  '/dist/bundle.css'
];

// Install event - cache static assets
self.addEventListener('install', (event) => {
  console.log('[SW] Installing service worker');
  event.waitUntil(
    caches.open(STATIC_CACHE).then((cache) => {
      return cache.addAll(STATIC_FILES);
    })
  );
  // Force activation of new service worker
  self.skipWaiting();
});

// Activate event - clean up old caches
self.addEventListener('activate', (event) => {
  console.log('[SW] Activating service worker');
  event.waitUntil(
    caches.keys().then((cacheNames) => {
      return Promise.all(
        cacheNames.map((cacheName) => {
          if (cacheName !== STATIC_CACHE && cacheName !== CACHE_NAME) {
            console.log('[SW] Deleting old cache:', cacheName);
            return caches.delete(cacheName);
          }
        })
      );
    })
  );
  // Take control of all clients immediately
  self.clients.claim();
});

// Fetch event - serve from cache, fallback to network
self.addEventListener('fetch', (event) => {
  const url = new URL(event.request.url);

  // Skip cross-origin requests
  if (url.origin !== location.origin) {
    return;
  }

  // Handle API requests - network first
  if (url.pathname.startsWith('/api/') ||
      url.pathname.startsWith('/object-molecules') ||
      url.pathname.startsWith('/generate-sdfs')) {
    event.respondWith(
      fetch(event.request)
        .then((response) => {
          // Cache successful API responses
          if (response.status === 200) {
            const responseClone = response.clone();
            caches.open(CACHE_NAME).then((cache) => {
              cache.put(event.request, responseClone);
            });
          }
          return response;
        })
        .catch(() => {
          // Fallback to cache for API requests
          return caches.match(event.request);
        })
    );
    return;
  }

  // Handle static assets and HTML - cache first, then network
  event.respondWith(
    caches.match(event.request)
      .then((response) => {
        if (response) {
          return response;
        }

        return fetch(event.request).then((response) => {
          // Don't cache non-successful responses
          if (!response || response.status !== 200 || response.type !== 'basic') {
            return response;
          }

          // Cache the response
          const responseClone = response.clone();
          caches.open(CACHE_NAME).then((cache) => {
            cache.put(event.request, responseClone);
          });

          return response;
        });
      })
  );
});

// Message event - handle cache invalidation from app
self.addEventListener('message', (event) => {
  if (event.data && event.data.type === 'SKIP_WAITING') {
    self.skipWaiting();
    // Send response to prevent runtime.lastError
    event.ports[0]?.postMessage({ success: true });
    return;
  }

  if (event.data && event.data.type === 'CLEAR_CACHE') {
    caches.keys().then((cacheNames) => {
      return Promise.all(
        cacheNames.map((cacheName) => {
          console.log('[SW] Clearing cache:', cacheName);
          return caches.delete(cacheName);
        })
      );
    }).then(() => {
      // Send response to prevent runtime.lastError
      event.ports[0]?.postMessage({ success: true });
    }).catch((error) => {
      // Send error response
      event.ports[0]?.postMessage({ success: false, error: error.message });
    });
    return;
  }

  // Send default response for any other messages
  event.ports[0]?.postMessage({ success: true, type: 'unknown' });
});
