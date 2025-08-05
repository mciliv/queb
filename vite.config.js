import { defineConfig } from 'vite'
import { resolve } from 'path'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  root: 'frontend',
  publicDir: 'assets',
  server: {
    port: 3000,
    proxy: {
      // Proxy API requests to our local Node.js backend
      '/api': {
        target: 'https://localhost:3002',
        changeOrigin: true,
        secure: false,
        rewrite: (path) => path.replace(/^\/api/, '')
      },
      // Proxy payment routes (not using /api prefix)
      '/stripe-config': 'https://localhost:3002',
      '/setup-payment-method': 'https://localhost:3002',
      '/update-payment-method': 'https://localhost:3002',
      // Proxy static file serving
      '/sdf_files': 'https://localhost:3002'
    }
  },
  build: {
    outDir: '../dist/frontend',
    emptyOutDir: true,
    rollupOptions: {
      input: {
        main: resolve(__dirname, 'frontend/index.html')
      }
    }
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, 'frontend')
    }
  }
}) 