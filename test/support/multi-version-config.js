// Multi-version testing configuration
// Enables running multiple versions/branches simultaneously on different ports

const VERSION_PORTS = {
  main: {
    http: 3000,
    https: 3001,
    debug: 9229,
    chrome: 9225
  },
  previous: {
    http: 3010, 
    https: 3011,
    debug: 9230,
    chrome: 9226
  },
  staging: {
    http: 3020,
    https: 3021, 
    debug: 9231,
    chrome: 9227
  },
  test: {
    http: 3333,
    https: 3334,
    debug: 9232,
    chrome: 9228
  }
};

const getVersionConfig = (version = 'main') => {
  const config = VERSION_PORTS[version];
  if (!config) {
    throw new Error(`Unknown version: ${version}. Available: ${Object.keys(VERSION_PORTS).join(', ')}`);
  }
  
  return {
    ...config,
    version,
    baseUrl: `http://localhost:${config.http}`,
    httpsUrl: `https://localhost:${config.https}`,
    debugUrl: `http://127.0.0.1:${config.chrome}`
  };
};

const getAllVersions = () => Object.keys(VERSION_PORTS);

const isPortAvailable = async (port) => {
  const net = require('net');
  return new Promise((resolve) => {
    const server = net.createServer();
    server.listen(port, (err) => {
      if (err) {
        resolve(false);
      } else {
        server.once('close', () => resolve(true));
        server.close();
      }
    });
    server.on('error', () => resolve(false));
  });
};

const checkVersionAvailability = async (version) => {
  const config = getVersionConfig(version);
  const ports = [config.http, config.https, config.debug, config.chrome];
  
  const availability = {};
  for (const port of ports) {
    availability[port] = await isPortAvailable(port);
  }
  
  return {
    version,
    available: Object.values(availability).every(Boolean),
    portStatus: availability,
    config
  };
};

module.exports = {
  VERSION_PORTS,
  getVersionConfig,
  getAllVersions,
  isPortAvailable,
  checkVersionAvailability
};
