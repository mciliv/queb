# Getting Started

## Prerequisites

- Node.js v18+
- Python 3.8+ with RDKit (`pip install rdkit-pypi`)
- OpenAI API key in `.env` as `OPENAI_API_KEY=sk-...`

Optional: PostgreSQL for user features, Google Cloud SDK for deployment.

## Setup

```
npm install
cp env.example .env   # then add your OPENAI_API_KEY
npm start             # runs at http://localhost:8080
```

## Key npm Scripts

| Script | Purpose |
|--------|---------|
| `npm start` | Start server |
| `npm run debug` | Start with debugger |
| `npm run build` | Build frontend |
| `npm test` | Run unit tests |
| `npm run setup:db` | Initialize database |
| `npm run cleanup:logs` | Clean old logs |

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `OPENAI_API_KEY` | Yes | OpenAI API key |
| `NODE_ENV` | No | `development` or `production` |
| `PORT` | No | Server port (default: 8080) |
| `DATABASE_URL` | No | PostgreSQL connection string |

## Verification

1. Type "water" → see H₂O in 3D
2. Camera mode → click object → see chemical analysis

## Troubleshooting

- **Module errors**: `rm -rf node_modules && npm install`
- **RDKit issues**: Use conda: `conda install -c conda-forge rdkit`
- **Port in use**: Change `PORT` in `.env`
