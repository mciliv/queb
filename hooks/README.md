# Git Hooks

This directory contains custom git hooks for the project.

## Pre-commit Hook

The `pre-commit` hook runs smoke tests before allowing commits.

### Behavior

- Runs `npm run test:smoke` to execute smoke tests
- Times out after 10 seconds if tests hang
- Allows commits to proceed even if tests fail or are unavailable
- This prevents commit stalling due to test environment issues

### Skipping Tests

If you need to skip the pre-commit tests temporarily:

```bash
git commit --no-verify
```

### Setup

The hooks path is configured automatically by the `prepare` script in package.json:

```json
{
  "scripts": {
    "prepare": "git config core.hooksPath hooks"
  }
}
```

Make sure the hooks are executable:
```bash
chmod +x hooks/pre-commit
```