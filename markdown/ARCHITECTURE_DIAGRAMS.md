# Architecture Diagrams - Molecular Space Analyzer

## 1. System Architecture Overview

```mermaid
graph TB
    subgraph "Client Layer"
        Browser[Web Browser]
        PWA[Progressive Web App]
        Camera[Device Camera]
    end
    
    subgraph "Frontend - React App"
        App[App.jsx<br/>Main Orchestrator]
        AppController[AppController.js<br/>State Management]
        
        subgraph "Input Components"
            TextInput[TextInput]
            CameraSection[CameraSection]
            PhotoSection[PhotoSection]
            LinkSection[LinkSection]
        end
        
        subgraph "Display Components"
            MolecularColumn[MolecularColumn]
            MoleculeViewer[MoleculeViewer<br/>3Dmol.js]
        end
        
        useApi[useApi Hook<br/>API Client]
    end
    
    subgraph "Backend - Express Server"
        Server[server.js<br/>Express API]
        
        subgraph "Services"
            Structuralizer[Structuralizer<br/>AI Analysis]
            MolecularProcessor[MolecularProcessor<br/>SDF Generation]
            NameResolver[NameResolver<br/>PubChem Client]
            UserService[UserService<br/>DB Management]
        end
        
        subgraph "Core Modules"
            Configuration[Configuration<br/>Config Manager]
            PromptEngine[PromptEngine<br/>Prompt Templates]
            ErrorHandler[ErrorHandler<br/>Error Management]
        end
    end
    
    subgraph "External Services"
        OpenAI[OpenAI API<br/>GPT-4o]
        PubChem[PubChem API<br/>Chemical Database]
    end
    
    subgraph "Data Storage"
        PostgreSQL[(PostgreSQL<br/>User Data)]
        SDF[SDF Files<br/>3D Structures]
        Logs[Log Files]
    end
    
    subgraph "Python Services"
        RDKit[RDKit<br/>sdf.py]
    end
    
    Browser --> App
    Camera --> CameraSection
    
    App --> AppController
    AppController --> TextInput
    AppController --> CameraSection
    AppController --> PhotoSection
    AppController --> LinkSection
    AppController --> MolecularColumn
    
    MolecularColumn --> MoleculeViewer
    
    TextInput --> useApi
    CameraSection --> useApi
    PhotoSection --> useApi
    LinkSection --> useApi
    
    useApi --> Server
    
    Server --> Structuralizer
    Server --> MolecularProcessor
    Server --> NameResolver
    Server --> UserService
    
    Structuralizer --> PromptEngine
    Structuralizer --> ErrorHandler
    Structuralizer --> Configuration
    Structuralizer --> OpenAI
    
    MolecularProcessor --> RDKit
    MolecularProcessor --> NameResolver
    MolecularProcessor --> SDF
    
    NameResolver --> PubChem
    
    UserService --> PostgreSQL
    
    Server --> Logs
    
    MoleculeViewer -.loads.-> SDF
    
    style OpenAI fill:#e1f5ff
    style PubChem fill:#e1f5ff
    style PostgreSQL fill:#ffe1e1
    style SDF fill:#ffe1e1
    style RDKit fill:#fff4e1
```

## 2. Data Flow - Text Analysis

```mermaid
sequenceDiagram
    participant User
    participant TextInput
    participant AppController
    participant useApi
    participant Server
    participant Structuralizer
    participant PromptEngine
    participant OpenAI
    participant NameResolver
    participant PubChem
    participant MolecularProcessor
    participant RDKit
    participant MoleculeViewer
    
    User->>TextInput: Enter "coffee"
    TextInput->>AppController: handleTextInput("coffee")
    AppController->>useApi: structuralizeText("coffee", "database")
    useApi->>Server: POST /structuralize
    
    Server->>Structuralizer: structuralizeText("coffee")
    Structuralizer->>PromptEngine: generateChemicalPrompt("coffee")
    PromptEngine-->>Structuralizer: Formatted prompt
    
    Structuralizer->>OpenAI: Chat completion (JSON mode)
    OpenAI-->>Structuralizer: {chemicals: [{name:"Water",smiles:null}, {name:"Caffeine",smiles:null}]}
    
    Note over Structuralizer: Enrich SMILES (database mode)
    
    loop For each chemical
        Structuralizer->>NameResolver: resolveName("Water")
        NameResolver->>PubChem: GET /compound/name/Water/property/SMILES
        PubChem-->>NameResolver: {smiles: "O"}
        NameResolver-->>Structuralizer: {smiles: "O"}
    end
    
    Structuralizer-->>Server: {chemicals: [{name:"Water",smiles:"O"}, {name:"Caffeine",smiles:"CN1C..."}]}
    Server-->>useApi: Response
    
    useApi->>Server: POST /generate-sdfs
    Server->>MolecularProcessor: processSmiles(["O", "CN1C..."])
    
    loop For each SMILES
        MolecularProcessor->>RDKit: sdf.py "O"
        alt RDKit succeeds
            RDKit-->>MolecularProcessor: /sdf_files/O.sdf
        else RDKit fails
            MolecularProcessor->>PubChem: Download SDF
            PubChem-->>MolecularProcessor: SDF content
        end
    end
    
    MolecularProcessor-->>Server: {sdfPaths: ["/sdf_files/O.sdf", ...]}
    Server-->>useApi: SDF paths
    
    useApi-->>AppController: Analysis complete
    AppController->>MolecularColumn: Display results
    MolecularColumn->>MoleculeViewer: Render 3D (sdfPath)
    MoleculeViewer->>MoleculeViewer: Load 3Dmol.js
    MoleculeViewer-->>User: Interactive 3D molecule
```

## 3. Data Flow - Camera Analysis

```mermaid
sequenceDiagram
    participant User
    participant Camera
    participant CameraSection
    participant AppController
    participant useApi
    participant Server
    participant Structuralizer
    participant OpenAI
    participant Flow[Same as Text Flow]
    
    User->>Camera: Grant permission
    Camera-->>CameraSection: Video stream
    
    User->>CameraSection: Click at (x, y)
    
    Note over CameraSection: Calculate crop region<br/>25% of min dimension
    
    CameraSection->>CameraSection: Capture frame to canvas
    CameraSection->>CameraSection: Extract crop region
    CameraSection->>CameraSection: Convert to base64
    
    CameraSection->>AppController: handleImageAnalysis()
    AppController->>useApi: structuralizeImage(base64, x, y, cropX, cropY, cropSize)
    
    useApi->>Server: POST /structuralize
    Note over Server: {imageBase64, x, y,<br/>cropMiddleX, cropMiddleY,<br/>cropSize}
    
    Server->>Structuralizer: structuralizeImage()
    Structuralizer->>OpenAI: Vision API + coordinates
    OpenAI-->>Structuralizer: {object: "coffee cup",<br/>chemicals: [...]}
    
    Note over Structuralizer,Flow: Continue with<br/>SMILES enrichment<br/>and SDF generation<br/>(same as text flow)
    
    Structuralizer-->>Server: Complete analysis
    Server-->>useApi: Response
    useApi-->>AppController: Results
    AppController->>User: Display molecules
```

## 4. Component Architecture

```mermaid
graph TB
    subgraph "Presentation Layer"
        UI[UI Components]
        State[State Management]
    end
    
    subgraph "Business Logic Layer"
        subgraph "Core Modules (Shared)"
            Config[Configuration]
            Prompt[PromptEngine]
            Error[ErrorHandler]
            Molecular[MolecularAnalysisService]
        end
    end
    
    subgraph "Service Layer"
        subgraph "Backend Services"
            Struct[Structuralizer]
            Processor[MolecularProcessor]
            Resolver[NameResolver]
            User[UserService]
        end
    end
    
    subgraph "Data Access Layer"
        API[API Client<br/>useApi]
        DB[Database<br/>PostgreSQL]
        Files[File System<br/>SDF Storage]
    end
    
    subgraph "External Integration Layer"
        AI[OpenAI Client]
        Chem[PubChem Client]
    end
    
    UI --> State
    State --> API
    
    API --> Struct
    API --> Processor
    
    Struct --> Config
    Struct --> Prompt
    Struct --> Error
    Struct --> Molecular
    Struct --> AI
    
    Processor --> Resolver
    Processor --> Files
    
    Resolver --> Chem
    
    User --> DB
    
    Config -.used by.-> Struct
    Config -.used by.-> Processor
    Config -.used by.-> Resolver
    
    style Config fill:#d4edda
    style Prompt fill:#d4edda
    style Error fill:#d4edda
    style Molecular fill:#d4edda
```

## 5. Class Diagram - Core Services

```mermaid
classDiagram
    class Configuration {
        -Map config
        +get(key) any
        +getDatabaseConfig() Object
        +validate() void
        +isDevelopment() boolean
        +isProduction() boolean
    }
    
    class PromptEngine {
        -Map templates
        -Map contexts
        +generateChemicalPrompt(object, options) string
        +generateDetectionPrompt(coords) string
        +generateNamePrompt(description) string
        +validateResponse(type, response) boolean
        +repairJSON(jsonString) Object
    }
    
    class ErrorHandler {
        -logger Logger
        +handle(error, context) HandledError
        +handleAIError(error, context) HandledError
        +handleValidationError(error, context) HandledError
        +initialize(logger) void
    }
    
    class Structuralizer {
        -client OpenAIClient
        -molecularProcessor MolecularProcessor
        +structuralize(payload) Object
        +structuralizeText(object, lookupMode) Object
        +structuralizeImage(imageData, coords) Object
    }
    
    class MolecularProcessor {
        -sdfDir string
        +processSmiles(smilesArray, overwrite) Object
        +generateSDF(smiles, overwrite) string
        +generateSDFByName(name, overwrite) Object
        -generateSmilesSDF(smiles) string
        -findExistingSdfFile(smiles) string
    }
    
    class NameResolver {
        -Map cache
        +resolveName(name) Object
        +downloadSDFBySmiles(smiles) string
    }
    
    class UserService {
        -pool PostgresPool
        +createUser(userData) Object
        +getUserByDeviceToken(token) Object
        +updateUser(token, data) void
        +incrementUsage(token) number
    }
    
    class useApi {
        -Map requestCache
        -Map pendingRequests
        +structuralizeText(text, mode) Promise
        +structuralizeImage(imageData, coords) Promise
        +generateSDFs(smiles, overwrite) Promise
        +nameToSdf(name, overwrite) Promise
    }
    
    Structuralizer --> PromptEngine : uses
    Structuralizer --> ErrorHandler : uses
    Structuralizer --> Configuration : uses
    Structuralizer --> MolecularProcessor : uses
    
    MolecularProcessor --> NameResolver : uses
    
    useApi --> Structuralizer : calls via HTTP
    useApi --> MolecularProcessor : calls via HTTP
```

## 6. State Machine - Input Mode Flow

```mermaid
stateDiagram-v2
    [*] --> Idle
    
    Idle --> TextMode: User types text
    Idle --> CameraMode: Click camera button
    Idle --> PhotoMode: Click photo button
    Idle --> LinkMode: Click link button
    
    TextMode --> Processing: Submit text
    CameraMode --> Processing: Click on video
    PhotoMode --> Processing: Upload file
    LinkMode --> Processing: Submit URL
    
    Processing --> AIAnalysis: Send to backend
    AIAnalysis --> SMILESEnrichment: Got chemical names
    
    SMILESEnrichment --> DatabaseLookup: mode=database
    SMILESEnrichment --> AILookup: mode=ai
    
    DatabaseLookup --> PubChemQuery
    AILookup --> OpenAIQuery
    
    PubChemQuery --> SDFGeneration: Got SMILES
    OpenAIQuery --> SDFGeneration: Got SMILES
    
    SDFGeneration --> RDKitPython: Try Python first
    RDKitPython --> DisplayResults: Success
    RDKitPython --> PubChemDownload: Fallback
    PubChemDownload --> DisplayResults: Got SDF
    
    DisplayResults --> Idle: Replace mode
    DisplayResults --> Accumulate: Accumulate mode
    Accumulate --> Idle: Continue analyzing
    
    Processing --> Error: Analysis failed
    Error --> Idle: User clears error
```

## 7. Deployment Architecture

```mermaid
graph TB
    subgraph "Client Devices"
        Desktop[Desktop Browser]
        Mobile[Mobile Browser]
        PWA[PWA Install]
    end
    
    subgraph "CDN Layer"
        CloudCDN[Google Cloud CDN]
    end
    
    subgraph "Google Cloud Platform"
        subgraph "App Engine"
            AppEngine1[Instance 1]
            AppEngine2[Instance 2]
            AppEngine3[Instance N]
        end
        
        subgraph "Cloud Functions"
            Function[molecular-analysis]
        end
        
        subgraph "Cloud Storage"
            StaticFiles[Static Assets<br/>SDF Files]
        end
        
        subgraph "Secret Manager"
            Secrets[API Keys<br/>Credentials]
        end
        
        subgraph "Cloud SQL"
            CloudSQL[(PostgreSQL<br/>User Data)]
        end
        
        LB[Load Balancer]
    end
    
    subgraph "External APIs"
        OpenAIAPI[OpenAI API]
        PubChemAPI[PubChem API]
    end
    
    Desktop --> CloudCDN
    Mobile --> CloudCDN
    PWA --> CloudCDN
    
    CloudCDN --> LB
    
    LB --> AppEngine1
    LB --> AppEngine2
    LB --> AppEngine3
    
    AppEngine1 --> StaticFiles
    AppEngine2 --> StaticFiles
    AppEngine3 --> StaticFiles
    
    AppEngine1 --> Secrets
    AppEngine1 --> CloudSQL
    AppEngine1 --> OpenAIAPI
    AppEngine1 --> PubChemAPI
    
    Function --> Secrets
    Function --> OpenAIAPI
    Function --> PubChemAPI
    
    style Desktop fill:#e3f2fd
    style Mobile fill:#e3f2fd
    style PWA fill:#e3f2fd
    style OpenAIAPI fill:#fff3e0
    style PubChemAPI fill:#fff3e0
```

## 8. Caching Architecture

```mermaid
graph LR
    subgraph "Frontend Caching"
        RequestCache[Request Cache<br/>5-30 min TTL]
        PendingRequests[Pending Requests<br/>Deduplication]
        BrowserCache[Browser Cache]
    end
    
    subgraph "Backend Caching"
        PubChemCache[PubChem Cache<br/>1 hour TTL<br/>LRU eviction]
        SDFFiles[SDF File Cache<br/>Permanent]
        PerformanceCache[Performance Cache<br/>In-Memory]
    end
    
    subgraph "Database Caching"
        MCPCache[MCP Cache Tables<br/>PostgreSQL]
    end
    
    Request[User Request] --> RequestCache
    RequestCache -->|Cache Miss| API[API Call]
    RequestCache -->|Cache Hit| Response[Return Cached]
    
    API --> PendingRequests
    PendingRequests -->|Identical Request| Wait[Wait for Pending]
    PendingRequests -->|New Request| Backend[Backend]
    
    Backend --> PubChemCache
    PubChemCache -->|Cache Miss| PubChemAPI[PubChem API]
    PubChemCache -->|Cache Hit| BackendResponse[Return]
    
    Backend --> SDFFiles
    SDFFiles -->|File Exists| BackendResponse
    SDFFiles -->|Not Found| Generate[Generate SDF]
    
    Generate --> SDFFiles
    
    MCPCache -.Optional.-> Backend
    
    style RequestCache fill:#c8e6c9
    style PubChemCache fill:#c8e6c9
    style SDFFiles fill:#c8e6c9
    style PerformanceCache fill:#c8e6c9
    style MCPCache fill:#c8e6c9
```

## 9. Error Handling Flow

```mermaid
graph TD
    Error[Error Occurs] --> Classify{Classify Error Type}
    
    Classify -->|AI Error| AIHandler[ErrorHandler.handleAIError]
    Classify -->|Validation Error| ValidationHandler[ErrorHandler.handleValidationError]
    Classify -->|Network Error| NetworkHandler[Network Error Handler]
    Classify -->|Unknown| GenericHandler[ErrorHandler.handle]
    
    AIHandler --> CheckRetry{Retryable?}
    NetworkHandler --> CheckRetry
    
    CheckRetry -->|Yes| RetryLogic[Retry with Backoff]
    CheckRetry -->|No| LogError[Log Error]
    
    RetryLogic --> AttemptCount{Attempts < Max?}
    AttemptCount -->|Yes| Retry[Retry Request]
    AttemptCount -->|No| LogError
    
    ValidationHandler --> LogError
    GenericHandler --> LogError
    
    LogError --> UserMessage[Generate User-Friendly Message]
    UserMessage --> Frontend[Send to Frontend]
    
    Frontend --> Display[Display Error UI]
    
    LogError -.development.-> Screenshot[Capture Screenshot]
    LogError --> FileLog[Write to Log File]
    
    style AIHandler fill:#ffccbc
    style ValidationHandler fill:#ffccbc
    style NetworkHandler fill:#ffccbc
    style GenericHandler fill:#ffccbc
    style LogError fill:#ef5350
    style Display fill:#90caf9
```

## 10. Module Dependency Graph

```mermaid
graph TD
    subgraph "Level 0 - No Dependencies"
        Logger[file-logger.js]
        Constants[constants.js]
        Patterns[validation-patterns.js]
    end
    
    subgraph "Level 1 - Core Modules"
        Config[Configuration.js]
        Prompt[PromptEngine.js]
        Error[ErrorHandler.js]
    end
    
    subgraph "Level 2 - Utility Services"
        NameRes[name-resolver.js]
        SecretMgr[secret-manager.js]
    end
    
    subgraph "Level 3 - Processing Services"
        MolProc[molecular-processor.js]
        UserSvc[user-service.js]
    end
    
    subgraph "Level 4 - Orchestration"
        Struct[Structuralizer.js]
        MolAnalysis[MolecularAnalysisService.js]
    end
    
    subgraph "Level 5 - API Layer"
        Server[server.js]
    end
    
    subgraph "Level 6 - Frontend"
        UseApi[useApi.js]
        AppCtrl[AppController.js]
        App[App.jsx]
    end
    
    Logger --> Error
    Logger --> Config
    
    Config --> Error
    Config --> Struct
    Config --> MolProc
    Config --> Server
    
    Prompt --> Struct
    Error --> Struct
    Error --> Server
    
    NameRes --> MolProc
    NameRes --> Struct
    
    MolProc --> Struct
    MolProc --> Server
    
    UserSvc --> Server
    
    Struct --> MolAnalysis
    Struct --> Server
    
    MolAnalysis --> Server
    
    Server --> UseApi
    UseApi --> AppCtrl
    AppCtrl --> App
    
    style Logger fill:#e8f5e9
    style Constants fill:#e8f5e9
    style Patterns fill:#e8f5e9
    style Config fill:#fff9c4
    style Prompt fill:#fff9c4
    style Error fill:#fff9c4
    style Server fill:#ffcdd2
    style App fill:#c5cae9
```

## How to View These Diagrams

These diagrams use **Mermaid** syntax and can be viewed in:

1. **GitHub/GitLab**: Renders automatically in `.md` files
2. **VS Code**: Install "Markdown Preview Mermaid Support" extension
3. **Cursor**: Should render in preview mode
4. **Online**: https://mermaid.live/
5. **Documentation sites**: Supports Mermaid (GitBook, Docusaurus, etc.)

## Diagram Export

Use `@mermaid-js/mermaid-cli` to export diagrams as PNG/SVG.

