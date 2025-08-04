# Vite + Node.js Development Setup

## Overview

We're using Vite for frontend development with Hot Module Replacement (HMR) and integrated Node.js backend with automatic reloading.

## How It Works

### Frontend (Vite)
- **Vite Dev Server**: Runs on port 3000 with HMR
- **Hot Module Replacement**: Instant updates without full page reload
- **Proxy Configuration**: API requests automatically proxied to backend
- **Fast Development**: Lightning-fast hot reloads

### Backend (Node.js)
- **Nodemon**: Automatically restarts server when backend files change
- **Integrated**: Runs alongside Vite in the same dev script
- **API Proxy**: Vite proxies API calls to backend on port 3001

## Ports Used

- **3000**: Vite dev server (frontend + HMR)
- **3001**: Node.js backend (API)
- **3002**: HTTPS server (production-like)

## Benefits

1. **Lightning Fast**: Vite's HMR is incredibly fast
2. **Integrated**: Both frontend and backend in one command
3. **Modern**: Uses latest frontend tooling
4. **Reliable**: No custom reload systems needed
5. **Proxy**: Seamless API communication between frontend and backend

## Usage

```bash
npm run dev
```

This starts:
- Vite dev server on port 3000 (frontend + HMR)
- Node.js backend on port 3001 (API)
- HTTPS server on port 3002

## Development URLs

- **Frontend**: `http://localhost:3000` (Vite + HMR)
- **Backend API**: `http://localhost:3001` (API endpoints)
- **HTTPS**: `https://dev.queb.space` (production-like)

## How to Use

1. **Frontend Changes**: Edit files in `frontend/` - changes appear instantly
2. **Backend Changes**: Edit files in `backend/` - server auto-restarts
3. **API Calls**: Frontend automatically proxies to backend

## Files

- `dev`: Main development script (Vite + Node.js)
- `vite.config.js`: Vite configuration with API proxy
- `cleanup`: Cleans up all processes 