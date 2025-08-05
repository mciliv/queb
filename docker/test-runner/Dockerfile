FROM node:20-alpine

# Install system dependencies
RUN apk add --no-cache \
    python3 \
    py3-pip \
    git \
    bash \
    curl \
    build-base \
    python3-dev

# Install Python dependencies
RUN pip3 install rdkit --break-system-packages

# Install PM2 globally
RUN npm install -g pm2

# Create app directory
WORKDIR /app

# Copy package files
COPY package*.json ./

# Install node dependencies
RUN npm ci --only=production

# Copy application files
COPY . .

# Create necessary directories
RUN mkdir -p logs sdf_files

# Create non-root user
RUN addgroup -g 1001 -S nodejs && \
    adduser -S nodejs -u 1001 && \
    chown -R nodejs:nodejs /app

USER nodejs

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD node scripts/health-check.js || exit 1

# Start PM2 with ecosystem config
CMD ["pm2-runtime", "start", "ecosystem.config.js"]