#!/usr/bin/env node

const path = require('path');
const fs = require('fs');
const os = require('os');

// Determine project directory
function findProjectRoot(startDir) {
  const scriptDir = __dirname;
  const projectRoot = path.dirname(scriptDir);
  
  if (fs.existsSync(path.join(projectRoot, 'package.json'))) {
    return projectRoot;
  }
  
  return projectRoot;
}

// Setup aliases function
function setupAliases(projectRoot) {
  const scriptPath = path.join(projectRoot, 'run');
  const homeDirectory = os.homedir();
  const aliasFilename = '.queb-aliases';
  const aliasFilePath = path.join(homeDirectory, aliasFilename);
  const aliasLines = [
    `alias r='${scriptPath}'`,
    `alias ship='${scriptPath} deploy'`,
    `alias clean='${scriptPath} clean'`,
    `alias test='${scriptPath} test'`,
    `alias build='${scriptPath} build'`,
    `alias dev='${scriptPath} dev'`,
    `alias start='${scriptPath} start'`
  ];
  const aliasContent = `${aliasLines.join('\n')}\n`;

  console.log(`🔧 Setting up development aliases for ${path.basename(projectRoot)}...`);

  const aliasFileExists = fs.existsSync(aliasFilePath);
  const aliasFileNeedsUpdate = !aliasFileExists || fs.readFileSync(aliasFilePath, 'utf8') !== aliasContent;

  if (aliasFileNeedsUpdate) {
    fs.writeFileSync(aliasFilePath, aliasContent, { mode: 0o644 });
    console.log(`✅ Aliases written to ${aliasFilePath}`);
  } else {
    console.log(`ℹ️ Aliases already up to date at ${aliasFilePath}`);
  }

  const shellCandidates = ['.zshrc', '.zprofile', '.bashrc', '.bash_profile', '.profile'];
  let shellConfig = shellCandidates.find((candidate) => fs.existsSync(path.join(homeDirectory, candidate)));
  if (!shellConfig) {
    shellConfig = '.zshrc';
  }
  const shellConfigPath = path.join(homeDirectory, shellConfig);
  const shellConfigContent = fs.existsSync(shellConfigPath) ? fs.readFileSync(shellConfigPath, 'utf8') : '';
  const sourceLine = `source "${aliasFilePath}"`;

  if (!shellConfigContent.includes(aliasFilePath)) {
    const needsNewline = shellConfigContent.length > 0 && !shellConfigContent.endsWith('\n');
    const prefix = needsNewline ? '\n' : '';
    fs.appendFileSync(shellConfigPath, `${prefix}${sourceLine}\n`);
    console.log(`✅ Added alias sourcing to ${shellConfigPath}`);
  } else {
    console.log(`ℹ️ ${shellConfigPath} already sources ${aliasFilePath}`);
  }

  console.log('\n🎉 Aliases configured. Restart your shell or run:');
  console.log(`  source ${shellConfigPath}`);
}

// Main execution
function main() {
  const projectRoot = findProjectRoot(process.cwd());
  setupAliases(projectRoot);
}

main();




