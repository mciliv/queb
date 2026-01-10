const { createContainer } = require('../../core/services');
const { startServer: startServerFromModule } = require('./server');

if (require.main === module) {
  const container = createContainer();
  startServerFromModule(container);
}

module.exports = {
  startServer: startServerFromModule,
};
