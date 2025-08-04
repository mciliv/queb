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
./run dev     # Start development server
./run server  # Start production server
./run test    # Run tests
./run format  # Format code
./run clean   # Clean processes
```

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