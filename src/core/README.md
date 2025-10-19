# Core Architecture Documentation

This directory contains the core modules that implement "A Philosophy of Software Design" principles throughout the molecular analysis application.

## Design Principles Applied

### 1. Deep Modules
Each core module provides a simple interface that hides complex implementation details:

- **Configuration**: Simple `get()` method hides complex environment variable loading and validation
- **PromptEngine**: Simple prompt generation methods hide complex template management and validation
- **ErrorHandler**: Simple `handle()` method hides complex error classification and recovery logic
- **MolecularAnalysisService**: Simple analysis methods hide complex AI integration and data processing

### 2. Information Hiding
Implementation details are carefully hidden from users:

- Configuration file formats and loading strategies
- Prompt template structures and validation rules
- Error classification algorithms and recovery strategies
- AI service integration complexities

### 3. General-Purpose Modules
Core modules are designed to be reusable across different parts of the application:

- Configuration can be used by any component needing environment variables
- ErrorHandler can handle errors from any source (AI, network, file system, etc.)
- PromptEngine can generate prompts for different AI tasks

## Module Overview

### Architecture Diagram

```mermaid
flowchart TB
  subgraph Core
    CFG[Configuration]\n- env loading\n- validation\n- type-safe getters
    ERR[ErrorHandler]\n- classify\n- message\n- recovery
    PENG[PromptEngine]\n- templates\n- validation\n- repair
    MAS[MolecularAnalysisService]\n- text/image analysis\n- enrichment\n- structures
  end

  %% Dependencies
  MAS --> PENG
  MAS --> ERR
  MAS --> CFG
  PENG --> ERR
  ERR --> CFG
```

### Configuration.js
**Purpose**: Unified configuration management with validation and type safety

**Key Features**:
- Hierarchical environment variable loading (.env.defaults → .env → .env.local)
- Type-safe getters with validation
- Environment-specific behavior (development vs production)
- Clear error messages for configuration issues

**Interface**:
```javascript
const config = require('./Configuration');

// Simple getters
const port = config.get('port');
const dbConfig = config.getDatabaseConfig();

// Environment checks
if (config.isProduction()) { /* ... */ }

// Validation
config.validate(); // Throws with clear error messages
```

### PromptEngine.js
**Purpose**: Centralized prompt management with validation and context awareness

**Key Features**:
- Template-based prompt generation
- Context-aware examples and formatting
- Response validation and repair
- Eliminates prompt duplication across the codebase

**Interface**:
```javascript
const promptEngine = require('./PromptEngine');

// Generate prompts
const prompt = promptEngine.generateChemicalPrompt('coffee');
const detection = promptEngine.generateDetectionPrompt({ x: 100, y: 200 });

// Validate responses
const isValid = promptEngine.validateResponse('chemical_analysis', response);

// Repair malformed JSON
const repaired = promptEngine.repairJSON(malformedJson);
```

### ErrorHandler.js
**Purpose**: Unified error handling with classification, logging, and recovery

**Key Features**:
- Automatic error classification and severity assessment
- User-friendly error message generation
- Recovery suggestion system
- Error frequency tracking for monitoring

**Interface**:
```javascript
const errorHandler = require('./ErrorHandler');

// Handle any error with context
const result = errorHandler.handle(error, { category: 'ai_service' });

// Specialized handlers
const aiResult = errorHandler.handleAIError(aiError);
const networkResult = errorHandler.handleNetworkError(networkError);

// Result contains user-friendly message and recovery info
console.log(result.message); // User-friendly message
console.log(result.recovery); // Recovery suggestion
```

### MolecularAnalysisService.js
**Purpose**: High-level service for molecular analysis tasks

**Key Features**:
- Unified interface for text and image analysis
- Automatic data enrichment and structure generation
- Dependency injection for testability
- Comprehensive error handling

**Interface**:
```javascript
const analysisService = require('./MolecularAnalysisService');

// Simple analysis methods
const textResult = await analysisService.analyzeText('coffee');
const imageResult = await analysisService.analyzeImage(imageData);
const structures = await analysisService.generateStructures(smilesArray);

// Service status
const status = analysisService.getStatus();
```

## Benefits of This Architecture

### Reduced Complexity
- Each module has a simple, focused interface
- Complex logic is hidden behind clean abstractions
- Dependencies are clearly defined and injected

### Improved Maintainability
- Changes to implementation don't affect interface users
- Error handling is consistent across the application
- Configuration is centralized and validated

### Better Testability
- Modules can be tested in isolation
- Dependencies can be mocked easily
- Error conditions can be simulated reliably

### Enhanced Reliability
- Comprehensive error handling with recovery
- Configuration validation prevents runtime issues
- Consistent logging and monitoring

## Usage Guidelines

### When to Use Core Modules

**Always use Configuration for**:
- Environment variables
- Application settings
- Feature flags
- Database connections

**Always use ErrorHandler for**:
- Exception handling
- User error messages
- Logging errors
- Recovery strategies

**Use PromptEngine when**:
- Generating AI prompts
- Validating AI responses
- Managing prompt templates

**Use MolecularAnalysisService for**:
- High-level analysis tasks
- Combining multiple services
- Business logic operations

### Best Practices

1. **Keep interfaces simple**: Don't expose implementation details
2. **Use dependency injection**: Make modules testable and flexible
3. **Handle errors consistently**: Always use ErrorHandler for exceptions
4. **Validate inputs**: Use Configuration for all environment variables
5. **Document intent**: Comments should explain why, not what

## Future Enhancements

- Add metrics collection to ErrorHandler
- Implement caching in PromptEngine
- Add configuration hot-reloading
- Extend MolecularAnalysisService with batch processing

This architecture provides a solid foundation that can evolve while maintaining clean interfaces and hiding complexity from users.


