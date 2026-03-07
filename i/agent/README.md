# UML Project Viewer

A web application that displays your project's code structure as a UML-style diagram with integrated `.cursorignore` management.

## Features

- **Visual Project Structure**: Displays all code files (`.js`, `.ts`, `.py`, `.jsx`, `.tsx`) with their classes and functions
- **Interactive Diagram**: Click any class or function to view its code
- **Cursorignore Management**: Toggle files in/out of `.cursorignore` with a simple checkbox
- **Real-time Updates**: Refresh button to reload the project structure

## Installation

```bash
npm install
```

## Usage

```bash
npm start
```

Then open your browser to `http://localhost:3000`

## How It Works

1. The server scans your project directory for code files
2. Parses JavaScript/TypeScript files using Acorn
3. Parses Python files with regex patterns
4. Displays the structure in an interactive web interface
5. Allows you to manage `.cursorignore` entries visually
6. Click any class/function to view the full source code

## Project Structure

- `server.js` - Express server with file parsing and API endpoints
- `public/index.html` - Main HTML interface
- `public/style.css` - Styling
- `public/app.js` - Frontend JavaScript for interactivity








