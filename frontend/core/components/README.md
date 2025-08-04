# Modular Component Architecture

This directory contains the modular components that break up the monolithic `index.html` into maintainable, focused pieces.

## Components

### `AppShell.js`
**Purpose**: Main application shell and keyboard shortcuts
**Responsibilities**:
- Keyboard shortcuts (Cmd+K/Ctrl+K for input focus)
- Text analysis handling
- Processing state management
- Payment section visibility
- Event coordination between components

**Key Features**:
- ✅ **Keyboard Shortcuts Preserved**: Cmd+K/Ctrl+K functionality maintained
- ✅ **Text Analysis**: Handles Enter key analysis
- ✅ **Processing States**: Shows/hides loading states
- ✅ **Error Handling**: Dispatches errors to ErrorHandler

### `InputSection.js`
**Purpose**: Handles all input methods (text, camera, photo)
**Responsibilities**:
- Text input management
- Camera mode switching
- Photo upload handling
- URL analysis
- Input state management

**Key Features**:
- ✅ **Mode Switching**: Camera/Photo mode toggle
- ✅ **File Upload**: Photo upload with drag & drop
- ✅ **URL Analysis**: Paste image URLs for analysis
- ✅ **Input States**: Enable/disable, focus, clear

### `MolecularViewer.js`
**Purpose**: 3D molecular visualization and SDF processing
**Responsibilities**:
- SDF file generation
- 3D molecular rendering
- Results display
- Chemical data extraction

**Key Features**:
- ✅ **3D Rendering**: Sphere representation with van der Waals radii
- ✅ **SDF Processing**: Converts SMILES to 3D structures
- ✅ **Results Display**: Grid layout with close buttons
- ✅ **Chemical Summary**: Statistics and skipped chemicals

### `ErrorHandler.js`
**Purpose**: Centralized error handling and display
**Responsibilities**:
- Error collection and display
- Error panel management
- Global error catching
- Error statistics

**Key Features**:
- ✅ **Error Panel**: Visual error display
- ✅ **Auto-cleanup**: Errors auto-hide after 10s
- ✅ **Global Catching**: Uncaught errors and promise rejections
- ✅ **Error Types**: Critical, warning, error classifications

## Architecture Benefits

### ✅ **Preserved Functionality**
- **Keyboard Shortcuts**: Cmd+K/Ctrl+K still works
- **All Input Methods**: Text, camera, photo upload
- **3D Visualization**: Molecular rendering intact
- **Error Handling**: Comprehensive error management

### ✅ **Modular Design**
- **Single Responsibility**: Each component has one clear purpose
- **Loose Coupling**: Components communicate via events
- **Easy Testing**: Components can be tested independently
- **Maintainable**: Changes isolated to specific components

### ✅ **Event-Driven Communication**
```javascript
// Components communicate via custom events
document.dispatchEvent(new CustomEvent('analysisResult', { detail: data }));
document.addEventListener('analysisResult', handler);
```

### ✅ **Clean Separation**
- **UI Logic**: Separated from business logic
- **State Management**: Each component manages its own state
- **Error Boundaries**: Errors contained within components
- **Reusability**: Components can be reused in different contexts

## Usage

### Main App Integration
```javascript
// app-modular.js coordinates all components
import { AppShell } from './components/AppShell.js';
import { MolecularViewer } from './components/MolecularViewer.js';
import { ErrorHandler } from './components/ErrorHandler.js';

class ModularMolecularApp {
  constructor() {
    this.appShell = new AppShell();
    this.molecularViewer = new MolecularViewer();
    this.errorHandler = new ErrorHandler();
  }
}
```

### Component Communication
```javascript
// Dispatch events for cross-component communication
document.dispatchEvent(new CustomEvent('appError', { 
  detail: { message: 'Error message', type: 'error' } 
}));

// Listen for events from other components
document.addEventListener('analysisResult', (e) => {
  const { output, objectName } = e.detail;
  // Handle analysis result
});
```

## Migration from Monolithic

The modular structure preserves all functionality from the original monolithic `index.html`:

1. **Keyboard Shortcuts**: ✅ Preserved in AppShell
2. **Text Analysis**: ✅ Preserved in AppShell + InputSection
3. **Camera Functionality**: ✅ Preserved via existing camera components
4. **Photo Upload**: ✅ Preserved in InputSection
5. **3D Visualization**: ✅ Preserved in MolecularViewer
6. **Error Handling**: ✅ Enhanced in ErrorHandler
7. **Payment System**: ✅ Preserved via existing payment components

## Development

### Adding New Components
1. Create component file in `components/` directory
2. Export class with clear responsibilities
3. Use event-driven communication
4. Add to `app-modular.js` initialization

### Testing Components
```javascript
// Test individual components
const inputSection = new InputSection();
inputSection.initialize();
// Test specific functionality
```

### Debugging
- Each component logs its initialization
- Error events provide detailed error information
- Components can be tested independently 