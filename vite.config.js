import { defineConfig } from 'vite'
import { resolve } from 'path'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  root: 'frontend',
  publicDir: 'assets',
  server: {
    port: 3000
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