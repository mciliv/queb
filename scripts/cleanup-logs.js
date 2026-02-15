#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

/**
 * Clean up log files at the start of development process
 * This runs before starting the server to ensure clean logs
 */

function cleanupLogs(projectRoot) {
  console.log('🧹 Cleaning up previous logs...');
  
  const logsDir = path.join(projectRoot, 'logs');
  let cleanedCount = 0;
  
  try {
    if (fs.existsSync(logsDir)) {
      const files = fs.readdirSync(logsDir);
      const logFiles = files.filter(file => file.endsWith('.log'));
      
      // Keep only the most recent 3 log files, delete the rest
      if (logFiles.length > 3) {
        const filesWithStats = logFiles.map(file => ({
          name: file,
          path: path.join(logsDir, file),
          mtime: fs.statSync(path.join(logsDir, file)).mtime
        }));
        
        // Sort by modification time (newest first)
        filesWithStats.sort((a, b) => b.mtime - a.mtime);
        
        // Delete files beyond the first 3 (oldest)
        const filesToDelete = filesWithStats.slice(3);
        
        filesToDelete.forEach(file => {
          try {
            fs.unlinkSync(file.path);
            cleanedCount++;
            console.log(`   🗑️  Removed old log: ${file.name}`);
          } catch (error) {
            console.log(`   ⚠️  Could not remove ${file.name}: ${error.message}`);
          }
        });
      }
      
      if (cleanedCount > 0) {
        console.log(`✅ Cleaned up ${cleanedCount} old log files`);
      } else {
        console.log('✅ Log directory is clean (≤3 files)');
      }
    } else {
      console.log('✅ No logs directory found');
    }
  } catch (error) {
    console.log(`⚠️  Log cleanup warning: ${error.message}`);
  }
  
  return cleanedCount;
}

// If run directly, clean logs in current project
if (require.main === module) {
  const projectRoot = process.cwd();
  cleanupLogs(projectRoot);
}

module.exports = { cleanupLogs };