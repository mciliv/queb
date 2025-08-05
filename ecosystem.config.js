module.exports = {
  apps: [
    {
      name: 'mol-test-watcher',
      script: './scripts/background-test-runner.js',
      args: 'watch',
      instances: 1,
      autorestart: true,
      watch: false,
      max_memory_restart: '500M',
      env: {
        NODE_ENV: 'development',
        OPENAI_API_KEY: process.env.OPENAI_API_KEY || 'test-api-key-for-testing'
      },
      error_file: './logs/test-watcher-error.log',
      out_file: './logs/test-watcher-out.log',
      log_file: './logs/test-watcher-combined.log',
      time: true,
      restart_delay: 5000,
      max_restarts: 10,
      min_uptime: '10s',
      cron_restart: '0 */6 * * *', // Restart every 6 hours to prevent memory leaks
    },
    {
      name: 'mol-dev-watcher',
      script: './scripts/dev-test-watcher.js',
      args: 'start',
      instances: 1,
      autorestart: true,
      watch: false,
      max_memory_restart: '300M',
      env: {
        NODE_ENV: 'development',
        OPENAI_API_KEY: process.env.OPENAI_API_KEY || 'test-api-key-for-testing'
      },
      error_file: './logs/dev-watcher-error.log',
      out_file: './logs/dev-watcher-out.log',
      log_file: './logs/dev-watcher-combined.log',
      time: true,
      restart_delay: 3000,
      max_restarts: 5,
      min_uptime: '5s',
    },
    {
      name: 'mol-test-scheduler',
      script: './scripts/test-scheduler.js',
      instances: 1,
      autorestart: true,
      watch: false,
      max_memory_restart: '200M',
      env: {
        NODE_ENV: 'production'
      },
      error_file: './logs/scheduler-error.log',
      out_file: './logs/scheduler-out.log',
      log_file: './logs/scheduler-combined.log',
      time: true,
      cron_restart: '0 0 * * *', // Daily restart at midnight
    }
  ],

  deploy: {
    production: {
      user: 'node',
      host: 'localhost',
      ref: 'origin/main',
      repo: 'git@github.com:mciliv/mol.git',
      path: '/var/www/mol',
      'post-deploy': 'npm install && pm2 reload ecosystem.config.js --env production'
    }
  }
};