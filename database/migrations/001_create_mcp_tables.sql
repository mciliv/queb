-- MCP Enhancement Tables
-- Run this to enable performance tracking and caching features

-- Performance metrics table
CREATE TABLE IF NOT EXISTS performance_metrics (
  id SERIAL PRIMARY KEY,
  operation VARCHAR(100) NOT NULL,
  duration_ms INTEGER NOT NULL,
  metadata JSONB DEFAULT '{}',
  timestamp TIMESTAMP DEFAULT NOW(),
  
  -- Indexes for common queries
  INDEX idx_performance_operation (operation),
  INDEX idx_performance_timestamp (timestamp),
  INDEX idx_performance_operation_time (operation, timestamp)
);

-- Analysis cache table
CREATE TABLE IF NOT EXISTS analysis_cache (
  id SERIAL PRIMARY KEY,
  operation VARCHAR(100) NOT NULL,
  cache_key VARCHAR(64) NOT NULL,
  input_text TEXT,
  result JSONB NOT NULL,
  hit_count INTEGER DEFAULT 0,
  created_at TIMESTAMP DEFAULT NOW(),
  
  -- Prevent duplicate cache entries
  UNIQUE(cache_key, operation),
  
  -- Indexes for lookups
  INDEX idx_cache_operation (operation),
  INDEX idx_cache_key (cache_key),
  INDEX idx_cache_created (created_at)
);

-- Validation failures tracking (optional)
CREATE TABLE IF NOT EXISTS validation_failures (
  id SERIAL PRIMARY KEY,
  pattern VARCHAR(200),
  input TEXT,
  timestamp TIMESTAMP DEFAULT NOW(),
  
  INDEX idx_validation_pattern (pattern),
  INDEX idx_validation_timestamp (timestamp)
);

-- SDF catalog (optional)
CREATE TABLE IF NOT EXISTS sdf_catalog (
  id SERIAL PRIMARY KEY,
  filename VARCHAR(255) UNIQUE NOT NULL,
  molecule_name VARCHAR(255),
  smiles TEXT,
  formula VARCHAR(100),
  mass DECIMAL(10, 4),
  file_path TEXT,
  file_size BIGINT,
  indexed_at TIMESTAMP DEFAULT NOW(),
  
  INDEX idx_sdf_molecule (molecule_name),
  INDEX idx_sdf_formula (formula)
);

-- Error log (optional, for enhanced error tracking)
CREATE TABLE IF NOT EXISTS error_log (
  id SERIAL PRIMARY KEY,
  error_id VARCHAR(100),
  error_code VARCHAR(50),
  message TEXT,
  severity VARCHAR(20),
  category VARCHAR(50),
  context JSONB,
  timestamp TIMESTAMP DEFAULT NOW(),
  
  INDEX idx_error_timestamp (timestamp),
  INDEX idx_error_severity (severity),
  INDEX idx_error_category (category)
);

-- Auto-cleanup old data (runs daily)
-- Performance metrics: keep 30 days
-- Cache: keep 7 days
-- Validation failures: keep 90 days
-- Error log: keep 30 days

CREATE OR REPLACE FUNCTION cleanup_old_metrics() 
RETURNS void AS $$
BEGIN
  DELETE FROM performance_metrics WHERE timestamp < NOW() - INTERVAL '30 days';
  DELETE FROM analysis_cache WHERE created_at < NOW() - INTERVAL '7 days';
  DELETE FROM validation_failures WHERE timestamp < NOW() - INTERVAL '90 days';
  DELETE FROM error_log WHERE timestamp < NOW() - INTERVAL '30 days';
END;
$$ LANGUAGE plpgsql;

-- Optional: Set up automatic cleanup (requires pg_cron extension)
-- SELECT cron.schedule('cleanup-metrics', '0 2 * * *', 'SELECT cleanup_old_metrics()');

COMMENT ON TABLE performance_metrics IS 'Tracks operation performance for monitoring and optimization';
COMMENT ON TABLE analysis_cache IS 'Caches expensive AI analysis results to save costs and improve speed';
COMMENT ON TABLE validation_failures IS 'Tracks validation pattern effectiveness for optimization';
COMMENT ON TABLE sdf_catalog IS 'Searchable catalog of molecular structure files';
COMMENT ON TABLE error_log IS 'Enhanced error tracking for monitoring and debugging';









