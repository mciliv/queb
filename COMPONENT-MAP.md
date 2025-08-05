# Quick Component Reference

## Camera & Image Processing
- **cameraManager** - Video feed, stream management, capture
- **cameraHandler** - Photo upload, base64 conversion, image processing

## Payment & Subscriptions  
- **paymentManager** - Full Stripe integration, subscription management
- **simplePaymentManager** - Simplified payment flow

## UI & Display
- **uiManager** - Result display, molecular rendering helpers, UI utilities
- **MolecularViewer** - 3D molecular visualization
- **InputSection** - All input methods (text, camera, photo)
- **ErrorHandler** - Error display and management
- **AppShell** - App shell, keyboard shortcuts

## Backend Services
- **AtomPredictor** - OpenAI vision API integration
- **MolecularProcessor** - SDF file generation, chemical processing
- **UserService** - User management, authentication
- **SimpleUserService** - Simplified user management

## DO NOT CREATE
❌ CameraModule (use cameraManager)
❌ ImageHandler (use cameraHandler)  
❌ PaymentProcessor (use paymentManager)
❌ UIHelper (use uiManager)
❌ VisionAPI (use AtomPredictor)
❌ ChemicalProcessor (use MolecularProcessor)