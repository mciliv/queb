#!/usr/bin/env node
"use strict";

try {
  const nodemon = require("nodemon");
  nodemon({
    script: "backend/api/server.js",
    watch: ["backend"],
    ext: "js json",
    ignore: ["frontend/**", "test/**"],
    exec: "node --inspect=9229",
    env: { ...process.env, PORT: "3000", NODE_ENV: "development" },
    stdout: true,
    quiet: true,
  });
  nodemon.on("quit", () => process.exit(0));
} catch (err) {
  console.error("nodemon not installed or failed to start (npm i -D nodemon)");
  process.exit(1);
}


