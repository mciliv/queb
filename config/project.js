module.exports = {
  name: 'molecular-analysis',
  description: 'Molecular analysis app with AI-powered chemical identification',
  domain: process.env.DOMAIN_NAME || '',
  dnsZone: process.env.DNS_ZONE_NAME || '',
  region: process.env.REGION || 'us-central1',
  functionName: process.env.FUNCTION_NAME || 'molecular-analysis',
  pythonVersion: process.env.PYTHON_VERSION || '3.12.9'
};
