#!/bin/bash

echo "🧪 Testing Component Connections"
echo "================================"

# Check if server is running
echo "🔍 Checking server status..."
if curl -s http://localhost:8080 > /dev/null; then
    echo "✅ Server is running on http://localhost:8080"
else
    echo "❌ Server not running. Starting server..."
    ./run server &
    sleep 3
fi

# Run file structure tests
echo ""
echo "📁 Testing file structure..."
node tests/test-runner.js

echo ""
echo "📝 Manual Testing Instructions:"
echo "1. Open http://localhost:8080 in your browser"
echo "2. Open browser console (F12)"
echo "3. Run: componentTests.runAllTests()"
echo "4. Run: cameraTests.runCameraTests()"
echo ""
echo "🎯 Specific Camera Test:"
echo "1. Click camera mode checkbox"
echo "2. Click on video feed"
echo "3. Grant camera permission"
echo "4. Click on object in camera view"
echo "5. Verify analysis results appear"
echo ""
echo "✅ Component connection testing complete!" 