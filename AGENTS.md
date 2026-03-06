# Agent Contributor Guide

This document outlines the engineering standards and workflows for AI agents (and human contributors) working on the 1 project.

## 1. Core Mandates

- **Security First:** Never log, print, or commit secrets, API keys, or sensitive credentials. Protect .env files.
- **Technical Integrity:** You are responsible for the entire lifecycle: implementation, testing, and validation.
- **Verification is Finality:** Never assume success. Run tests and project-specific build commands (e.g., npm run build) after changes.

## 2. Engineering Standards

- **Coding Style:** Adhere to existing naming, formatting, typing, and commenting conventions. Analyze surrounding code before modifying.
- **Code Reusability:** Favor composition over duplication. Extract shared logic into src/core or src/utils.
- **Naming over Comments:** Prioritize self-documenting code with clear variable and function names.
- **Shell Scripts:** All scripts must be compatible with both bash and zsh.

## 3. Architecture & DI

The project uses a custom dependency injection container in src/core/services.js.
- **Controllers:** Keep them thin (HTTP <=> Domain translation only).
- **Use Cases:** Abstract domain logic into use-cases in src/server/services/.
- **Adapters:** Wrap infrastructure in adapters to match domain contracts.

## 4. Testing Requirements

- **Unit Tests:** All core logic must have unit tests in tests/suites/unit/.
- **Integration Tests:** Use tests/suites/integration/ for end-to-end service testing.
- **Regression:** Reproduce bugs with a test case before applying the fix.
- **Validation:** Run npm run test before submitting changes.

## 5. Development Workflow

1.  Research: Map the codebase and validate assumptions with grep and find.
2.  Strategy: Formulate a plan before execution.
3.  Execute: Apply surgical, targeted changes.
4.  Validate: Run tests and build the application (npm run build).

## 6. Ubuntu VM

Source: `vm.swift`. Run `./ubuntu-vm -h` for full usage, build instructions, and configuration details.
