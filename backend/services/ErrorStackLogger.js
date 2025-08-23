class ErrorStackLogger {
  static log(error, label = '') {
    if (!error) return;
    const prefix = label ? `${label}: ` : '';
    if (error.stack) {
      console.error(prefix + error.stack);
    } else if (error.message) {
      console.error(prefix + error.message);
    } else {
      console.error(prefix + String(error));
    }
  }
}

module.exports = ErrorStackLogger;


