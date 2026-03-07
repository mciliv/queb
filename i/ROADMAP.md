# Roadmap

## Updater

Background auto-update service for the CLI agent.

### Responsibilities
- Check for updates on startup (non-blocking)
- Detect when the running process is idle or between tasks (no active tool execution, no pending AI response) to avoid conflicts
- Download and apply updates in the background
- On next prompt or startup, print a single line summarizing what changed (e.g. `updated: v0.2.0 — added search tool, fixed shell escaping`)

### Design constraints
- Never interrupt an active agent loop or tool execution
- No confirmation prompts — updates are silent and automatic
- Update check and download happen in a background worker/child process
- The updated code takes effect on the next CLI invocation, not mid-session
- Integrity: verify checksums or signatures before applying

### Implementation sketch
1. `Updater` class with `checkForUpdate()`, `download()`, `apply()` methods
2. On startup, spawn a detached background process that checks the registry/repo for a newer version
3. If a new version is available, download to a staging directory and verify integrity
4. Swap the staged files into place when no agent loop is running (between invocations or at clean exit)
5. Write a `.update-log` file with the changelog summary
6. On next startup, if `.update-log` exists, print its contents and delete the file
