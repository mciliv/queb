# Auto-Commit Message Setup

## Using AI-Generated Commit Messages (Optional)

To enable AI-generated commit messages using OpenAI:

1. Set your OpenAI API key:
   ```bash
   export OPENAI_API_KEY="your-api-key-here"
   ```

2. Or add it to your `.env` file:
   ```
   OPENAI_API_KEY=your-api-key-here
   ```

## Fallback Mode

If no OpenAI API key is configured, the script will automatically use a simple commit message generator that analyzes your changes and creates conventional commit messages based on:
- File types changed
- Number of files
- Directory structure
- Common patterns (test files, docs, etc.)

## Usage

Run the development server with auto-commit:
```bash
npm run dev:commit
```

When you stop the server (Ctrl+C), it will:
1. Show you the changes
2. Generate a commit message
3. Ask for confirmation (y/n/e to edit)
4. Commit and push if confirmed
