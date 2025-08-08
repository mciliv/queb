---
alwaysApply: true
---

# Claude-4 instructions
Never assume perfection
Maximally simple & straightforward
only respond in code
If anything can be seen as a principle, put here
No logging
Prefer code comments over chat text
Show code chunks not descriptions
Don't create extraneous files; no markdown, if anything, simplify the code & tell me why to look at it

After completing change, run `dev` in an independent process
If all tests pass,
```sh
git commit -am "Brief description"
git push origin main
```

# Conventions

Name files, variables, functions, classes specifically