# Dependabot Auto-Merge Setup

## 1. Enable Auto-Merge in Repository Settings

Go to your repository on GitHub:
1. Settings → General
2. Scroll to "Pull Requests"
3. Enable "Allow auto-merge"
4. Enable "Automatically delete head branches"

## 2. Set Branch Protection Rules

1. Settings → Branches
2. Add rule for `main` or `master`
3. Enable:
   - Require pull request reviews before merging
   - Dismiss stale pull request approvals
   - Require status checks to pass (select your CI workflow)
   - Include administrators

## 3. GitHub CLI Auto-Merge (Alternative)

For immediate auto-merge of current Dependabot PRs:

```bash
# List all Dependabot PRs
gh pr list --author dependabot

# Auto-merge a specific PR
gh pr merge [PR-NUMBER] --auto --squash

# Auto-merge all Dependabot security PRs
gh pr list --author dependabot --search "security in:title" --json number --jq '.[].number' | xargs -I {} gh pr merge {} --auto --squash
```

## 4. Manual Security Update

To manually update dependencies with security issues:

```bash
# Check for vulnerabilities
npm audit

# Auto-fix vulnerabilities
npm audit fix

# Force fixes (use with caution)
npm audit fix --force
```
