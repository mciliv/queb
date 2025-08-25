#!/bin/bash
# Test script - run all test suites
npx jest --config test/jest.config.js --testPathPattern=unit.test.js --verbose --silent && npx jest --config test/jest.config.js --testPathPattern=integration.test.js --verbose && npx jest --config test/jest.config.js --testPathPattern=system.test.js --verbose --detectOpenHandles 