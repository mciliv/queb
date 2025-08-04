import { defineConfig } from 'vite'
import { resolve } from 'path'

export default defineConfig({
  root: 'frontend',
  publicDir: 'assets',
  server: {
    port: 3000,
    proxy: {
      // Proxy API requests to our Node.js backend
      '/api': {
        target: 'http://localhost:3001',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      },
      // Proxy other backend routes
      '/stripe-config': 'http://localhost:3001',
      '/setup-payment-method': 'http://localhost:3001',
      '/update-payment-method': 'http://localhost:3001',
      '/analyze-image': 'http://localhost:3001',
      '/sdf_files': 'http://localhost:3001'
    }
  },
  build: {
    outDir: '../dist/frontend',
    emptyOutDir: true,
    rollupOptions: {
      input: {
        main: resolve(__dirname, 'frontend/core/index.html'),
        app: resolve(__dirname, 'frontend/core/app.js')
      }
    }
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, 'frontend')
    }
  }
}) 