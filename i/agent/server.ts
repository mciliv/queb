import express, { Request, Response } from 'express';
import fs from 'fs/promises';
import path from 'path';
import { fileURLToPath } from 'url';
import { parse } from 'acorn';
import * as walk from 'acorn-walk';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const app = express();
const PORT = 3000;

// Type definitions
interface ClassInfo {
  name: string;
  type: 'class' | 'function';
  line: number;
}

interface FileInfo {
  path: string;
  name: string;
  fullPath: string;
  classes: ClassInfo[];
}

interface CursorignoreRequest {
  paths: string[];
  action: 'add' | 'remove';
  dir?: string;
}

// Helper to resolve directory path (handles ~ expansion)
function resolveDir(dir: string): string {
  if (dir.startsWith('~')) {
    return path.join(process.env.HOME || '', dir.slice(1));
  }
  return path.resolve(dir);
}

app.use(express.json());
app.use(express.static('public'));

// Parse JavaScript/TypeScript files to extract classes
async function parseJSFile(filePath: string): Promise<ClassInfo[]> {
  const content = await fs.readFile(filePath, 'utf-8');
  const classes: ClassInfo[] = [];
  
  try {
    const ast = parse(content, { 
      ecmaVersion: 'latest' as any, 
      sourceType: 'module',
      locations: true 
    });
    
    walk.simple(ast, {
      ClassDeclaration(node: any) {
        classes.push({
          name: node.id.name,
          type: 'class',
          line: node.loc.start.line
        });
      },
      FunctionDeclaration(node: any) {
        if (node.id) {
          classes.push({
            name: node.id.name,
            type: 'function',
            line: node.loc.start.line
          });
        }
      }
    });
  } catch (e) {
    // If parsing fails, just return file info
    console.log(`Error parsing ${filePath}:`, (e as Error).message);
  }
  
  return classes;
}

// Parse Python files to extract classes
async function parsePythonFile(filePath: string): Promise<ClassInfo[]> {
  const content = await fs.readFile(filePath, 'utf-8');
  const classes: ClassInfo[] = [];
  const lines = content.split('\n');
  
  lines.forEach((line, index) => {
    const classMatch = line.match(/^class\s+(\w+)/);
    const funcMatch = line.match(/^def\s+(\w+)/);
    
    if (classMatch) {
      classes.push({
        name: classMatch[1],
        type: 'class',
        line: index + 1
      });
    } else if (funcMatch) {
      classes.push({
        name: funcMatch[1],
        type: 'function',
        line: index + 1
      });
    }
  });
  
  return classes;
}

// Recursively scan directory for code files
async function scanDirectory(dir: string, baseDir: string = dir): Promise<FileInfo[]> {
  const files: FileInfo[] = [];
  const items = await fs.readdir(dir, { withFileTypes: true });
  
  for (const item of items) {
    const fullPath = path.join(dir, item.name);
    const relativePath = path.relative(baseDir, fullPath);
    
    // Skip node_modules and hidden directories
    if (item.name.startsWith('.') || item.name === 'node_modules' || item.name === 'public') {
      continue;
    }
    
    if (item.isDirectory()) {
      const subFiles = await scanDirectory(fullPath, baseDir);
      files.push(...subFiles);
    } else if (item.isFile()) {
      const ext = path.extname(item.name);
      if (['.js', '.ts', '.py', '.jsx', '.tsx'].includes(ext)) {
        let classes: ClassInfo[] = [];
        
        if (['.js', '.ts', '.jsx', '.tsx'].includes(ext)) {
          classes = await parseJSFile(fullPath);
        } else if (ext === '.py') {
          classes = await parsePythonFile(fullPath);
        }
        
        files.push({
          path: relativePath,
          name: item.name,
          fullPath: fullPath,
          classes: classes
        });
      }
    }
  }
  
  return files;
}

// API endpoint to get project structure
app.get('/api/project-structure', async (req: Request, res: Response) => {
  try {
    const targetDir = resolveDir((req.query.dir as string) || __dirname);
    console.log(`Scanning directory: ${targetDir}`);
    const files = await scanDirectory(targetDir);
    res.json({ files });
  } catch (error) {
    res.status(500).json({ error: (error as Error).message });
  }
});

// API endpoint to get file content
app.get('/api/file-content', async (req: Request, res: Response) => {
  try {
    const { path: filePath, dir } = req.query as { path: string; dir?: string };
    const targetDir = resolveDir(dir || __dirname);
    const fullPath = path.join(targetDir, filePath);
    const content = await fs.readFile(fullPath, 'utf-8');
    res.json({ content });
  } catch (error) {
    res.status(500).json({ error: (error as Error).message });
  }
});

// API endpoint to update .cursorignore
app.post('/api/cursorignore', async (req: Request<{}, {}, CursorignoreRequest>, res: Response) => {
  try {
    const { paths, action, dir } = req.body;
    const targetDir = resolveDir(dir || __dirname);
    const cursorignorePath = path.join(targetDir, '.cursorignore');
    
    let content = '';
    try {
      content = await fs.readFile(cursorignorePath, 'utf-8');
    } catch (e) {
      // File doesn't exist yet
    }
    
    let lines = content.split('\n').filter(line => line.trim());
    
    if (action === 'add') {
      paths.forEach(p => {
        if (!lines.includes(p)) {
          lines.push(p);
        }
      });
    } else if (action === 'remove') {
      lines = lines.filter(line => !paths.includes(line));
    }
    
    await fs.writeFile(cursorignorePath, lines.join('\n') + '\n');
    res.json({ success: true });
  } catch (error) {
    res.status(500).json({ error: (error as Error).message });
  }
});

// API endpoint to get current .cursorignore contents
app.get('/api/cursorignore', async (req: Request, res: Response) => {
  try {
    const targetDir = resolveDir((req.query.dir as string) || __dirname);
    const cursorignorePath = path.join(targetDir, '.cursorignore');
    let content = '';
    try {
      content = await fs.readFile(cursorignorePath, 'utf-8');
    } catch (e) {
      // File doesn't exist yet
    }
    const paths = content.split('\n').filter(line => line.trim());
    res.json({ paths });
  } catch (error) {
    res.status(500).json({ error: (error as Error).message });
  }
});

app.listen(PORT, () => {
  console.log(`Server running at http://localhost:${PORT}`);
});
