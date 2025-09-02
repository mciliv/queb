// Minimal no-op logger to satisfy imports in production build
const noop = () => {};

const logger = {
  info: (...args) => noop(args),
  warn: (...args) => noop(args),
  error: (...args) => noop(args),
  debug: (...args) => noop(args)
};

export default logger;

