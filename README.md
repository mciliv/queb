# Mol - Molecular Visualization App

Simple molecular analysis and 3D visualization tool.

## Setup

1. Install dependencies:
```bash
npm install
```

2. Create `.env` file:
```bash
OPENAI_API_KEY=your_key_here
NODE_ENV=development
PORT=8080
```

3. Optional - PostgreSQL database:
```bash
createdb mol_users
```

## Usage

```bash
dev              # Start development server
server           # Start production server
test             # Run tests
test:connections # Test backwards endpoint connections
format           # Format code
clean            # Clean processes
```

## Connection Testing

The app includes backwards connection testing that validates the entire data flow:

**Test Flow (backwards from output):**
1. Frontend Display ← Can we serve the frontend?
2. SDF File Serving ← Can we serve generated SDF files?  
3. SDF Generation ← Can we generate SDF from SMILES?
4. Molecular Analysis ← Can we analyze text/images to get SMILES?
5. Input Validation ← Can we handle valid inputs?


- `./run test:connections` - Run backwards connection tests
- `/health/connections` - Real-time connection health API
- `Ctrl+Shift+H` - Toggle connection health display (dev mode)

## Features

- Text input molecular analysis
- Camera/photo molecular analysis
- 3D molecular visualization (spheres only, van der Waals radii)
- Simple, minimal UI
- Optional payment system

## Tech Stack

- Frontend: Vanilla JS, 3Dmol.js, Vite
- Backend: Node.js, Express
- AI: OpenAI Vision API
- Database: PostgreSQL (optional)