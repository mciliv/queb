const path = require('path');

class OneLineReporter {
  constructor(globalConfig, options) {
    this._globalConfig = globalConfig;
    this._options = options || {};
  }

  onRunStart() {
    // no-op
  }

  onTestCaseResult(test, testCaseResult) {
    const status = testCaseResult.status || 'unknown';
    const icon = status === 'passed' ? 'PASS' : status === 'failed' ? 'FAIL' : status === 'skipped' ? 'SKIP' : status === 'todo' ? 'TODO' : status.toUpperCase();

    const relFile = path.relative(process.cwd(), test.path || '');
    const suite = (testCaseResult.ancestorTitles || []).join(' > ');
    const title = testCaseResult.title || '';
    const ms = Number.isFinite(testCaseResult.duration) ? `${testCaseResult.duration}ms` : '-';

    const line = `${icon} ${relFile}${suite ? ' > ' + suite : ''} > ${title} (${ms})`;
    process.stdout.write(line + '\n');
  }

  onTestResult() {
    // no-op: per-test lines already printed in onTestCaseResult
  }

  onRunComplete() {
    // Keep Jest default summary separate; this reporter is purely per-test-line
  }
}

module.exports = OneLineReporter;
