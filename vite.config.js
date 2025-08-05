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
      // Proxy other backend routes to local server
      '/stripe-config': 'https://localhost:3002',
      '/setup-payment-method': 'https://localhost:3002',
      '/update-payment-method': 'https://localhost:3002',
      '/analyze-image': 'https://localhost:3002',
      '/analyze-text': 'https://localhost:3002',
      '/sdf_files': 'https://localhost:3002'
    }
  },
  build: {
    outDir: '../dist/frontend',
    emptyOutDir: true,
    rollupOptions: {
      input: {
        main: resolve(__dirname, 'frontend/core/index.html')
      }
    }
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, 'frontend')
    }
  }
}) 