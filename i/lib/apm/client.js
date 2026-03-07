const http = require('http');

const PORT = process.env.APM_PORT || 7243;
const HOST = '127.0.0.1';
let _id = null;

function getId() {
  if (!_id) _id = `agent-${Math.random().toString(36).substr(2, 9)}`;
  return _id;
}

function report(data) {
  const payload = JSON.stringify(data);
  const req = http.request({
    hostname: HOST, port: PORT,
    path: `/ingest/${getId()}`, method: 'POST',
    headers: { 'Content-Type': 'application/json', 'Content-Length': Buffer.byteLength(payload) }
  }, (res) => res.resume());
  req.on('error', () => {});
  req.end(payload);
}

function getDirectives() {
  return new Promise(resolve => {
    http.get(`http://${HOST}:${PORT}/directive/${getId()}`, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => { try { resolve(JSON.parse(d)); } catch { resolve([]); } });
    }).on('error', () => resolve([]));
  });
}

function sendDirective(agentId, message) {
  const payload = JSON.stringify({ message, from: getId() });
  const req = http.request({
    hostname: HOST, port: PORT,
    path: `/directive/${agentId}`, method: 'POST',
    headers: { 'Content-Type': 'application/json', 'Content-Length': Buffer.byteLength(payload) }
  }, (res) => res.resume());
  req.on('error', () => {});
  req.end(payload);
}

function fetchAgents() {
  return new Promise(resolve => {
    http.get(`http://${HOST}:${PORT}/agents`, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => { try { resolve(JSON.parse(d)); } catch { resolve([]); } });
    }).on('error', () => resolve([]));
  });
}

module.exports = { report, getId, getDirectives, sendDirective, fetchAgents };
