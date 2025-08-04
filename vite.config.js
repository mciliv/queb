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
        target: 'https://dev.queb.space',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      },
      // Proxy other backend routes
      '/stripe-config': 'https://dev.queb.space',
      '/setup-payment-method': 'https://dev.queb.space',
      '/update-payment-method': 'https://dev.queb.space',
      '/analyze-image': 'https://dev.queb.space',
      '/sdf_files': 'https://dev.queb.space'
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