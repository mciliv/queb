#!/bin/bash
# Test script - run all test suites
npx jest --testPathPattern=unit.test.js --verbose --silent && npx jest --testPathPattern=integration.test.js --verbose && npx jest --testPathPattern=system.test.js --verbose --detectOpenHandles 