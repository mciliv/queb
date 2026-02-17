// Minimal console logger for services.js startup messages
// Prompt: services.js line 4 imports { log, warn, error } for startup logging
function log(message, data) {
  if (data !== undefined) console.log(message, data);
  else console.log(message);
}

function warn(message, data) {
  if (data !== undefined) console.warn(message, data);
  else console.warn(message);
}

function error(message, data) {
  if (data !== undefined) console.error(message, data);
  else console.error(message);
}

module.exports = { log, warn, error };
