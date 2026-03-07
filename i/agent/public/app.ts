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

interface ProjectStructureResponse {
    files: FileInfo[];
    error?: string;
}

interface FileContentResponse {
    content: string;
    error?: string;
}

interface CursorignoreResponse {
    paths: string[];
    error?: string;
}

interface SuccessResponse {
    success: boolean;
    error?: string;
}

// Global state
let projectData: FileInfo[] | null = null;
let ignoredPaths: string[] = [];
let activeElement: HTMLElement | null = null;

async function loadProjectStructure(): Promise<void> {
    const directoryInput = document.getElementById('directory-input') as HTMLInputElement;
    const directory = directoryInput.value.trim();
    
    if (!directory) {
        alert('Please enter a directory path');
        return;
    }
    
    try {
        const [structureRes, ignoreRes] = await Promise.all([
            fetch(`/api/project-structure?dir=${encodeURIComponent(directory)}`),
            fetch(`/api/cursorignore?dir=${encodeURIComponent(directory)}`)
        ]);
        
        const structureData: ProjectStructureResponse = await structureRes.json();
        const ignoreData: CursorignoreResponse = await ignoreRes.json();
        
        if (structureData.error) {
            throw new Error(structureData.error);
        }
        
        projectData = structureData.files;
        ignoredPaths = ignoreData.paths;
        
        renderDiagram();
    } catch (error) {
        console.error('Error loading project structure:', error);
        const container = document.getElementById('diagram-container');
        if (container) {
            container.innerHTML = 
                `<div class="empty-state">Error: ${(error as Error).message}</div>`;
        }
    }
}

function renderDiagram(): void {
    const container = document.getElementById('diagram-container');
    if (!container) return;
    
    if (!projectData || projectData.length === 0) {
        container.innerHTML = '<div class="empty-state">No code files found in the project</div>';
        return;
    }
    
    container.innerHTML = '';
    
    projectData.forEach(file => {
        const fileGroup = document.createElement('div');
        fileGroup.className = 'file-group';
        
        if (ignoredPaths.includes(file.path)) {
            fileGroup.classList.add('ignored');
        }
        
        const fileHeader = document.createElement('div');
        fileHeader.className = 'file-header';
        
        const fileName = document.createElement('div');
        fileName.className = 'file-name';
        fileName.textContent = file.path;
        
        const cursorignoreToggle = document.createElement('label');
        cursorignoreToggle.className = 'cursorignore-toggle';
        
        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.checked = ignoredPaths.includes(file.path);
        checkbox.addEventListener('change', (e) => {
            const target = e.target as HTMLInputElement;
            toggleCursorignore(file.path, target.checked);
        });
        
        const label = document.createElement('span');
        label.textContent = 'Ignore in Cursor';
        
        cursorignoreToggle.appendChild(checkbox);
        cursorignoreToggle.appendChild(label);
        
        fileHeader.appendChild(fileName);
        fileHeader.appendChild(cursorignoreToggle);
        
        const classesContainer = document.createElement('div');
        classesContainer.className = 'classes-container';
        
        if (file.classes.length === 0) {
            const emptyMsg = document.createElement('div');
            emptyMsg.textContent = 'No classes or functions found';
            emptyMsg.style.color = '#858585';
            emptyMsg.style.fontSize = '0.85rem';
            classesContainer.appendChild(emptyMsg);
        } else {
            file.classes.forEach(cls => {
                const classBox = document.createElement('div');
                classBox.className = 'class-box';
                classBox.dataset.file = file.path;
                classBox.dataset.line = cls.line.toString();
                
                const className = document.createElement('div');
                className.className = 'class-name';
                className.textContent = cls.name;
                
                const classType = document.createElement('div');
                classType.className = 'class-type';
                classType.textContent = cls.type;
                
                classBox.appendChild(className);
                classBox.appendChild(classType);
                
                classBox.addEventListener('click', () => {
                    openFile(file.path, classBox);
                });
                
                classesContainer.appendChild(classBox);
            });
        }
        
        fileGroup.appendChild(fileHeader);
        fileGroup.appendChild(classesContainer);
        container.appendChild(fileGroup);
    });
}

async function toggleCursorignore(filePath: string, isIgnored: boolean): Promise<void> {
    const directoryInput = document.getElementById('directory-input') as HTMLInputElement;
    const directory = directoryInput.value.trim();
    
    try {
        const action = isIgnored ? 'add' : 'remove';
        await fetch('/api/cursorignore', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ paths: [filePath], action, dir: directory })
        });
        
        if (isIgnored) {
            ignoredPaths.push(filePath);
        } else {
            ignoredPaths = ignoredPaths.filter(p => p !== filePath);
        }
        
        renderDiagram();
    } catch (error) {
        console.error('Error updating .cursorignore:', error);
    }
}

async function openFile(filePath: string, element: HTMLElement): Promise<void> {
    const directoryInput = document.getElementById('directory-input') as HTMLInputElement;
    const directory = directoryInput.value.trim();
    
    try {
        // Update active element
        if (activeElement) {
            activeElement.classList.remove('active');
        }
        activeElement = element;
        element.classList.add('active');
        
        const res = await fetch(`/api/file-content?path=${encodeURIComponent(filePath)}&dir=${encodeURIComponent(directory)}`);
        const data: FileContentResponse = await res.json();
        
        const fileNameEl = document.getElementById('file-name');
        const codeContentEl = document.getElementById('code-content');
        
        if (fileNameEl) fileNameEl.textContent = filePath;
        if (codeContentEl) codeContentEl.textContent = data.content;
    } catch (error) {
        console.error('Error loading file:', error);
        const codeContentEl = document.getElementById('code-content');
        if (codeContentEl) codeContentEl.textContent = 'Error loading file content';
    }
}

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    const refreshBtn = document.getElementById('refresh-btn');
    if (refreshBtn) {
        refreshBtn.addEventListener('click', loadProjectStructure);
    }
    loadProjectStructure();
});
