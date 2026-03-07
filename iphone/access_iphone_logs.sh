#!/bin/bash

# Simple iPhone Log Access Script
# Run this when your iPhone is connected via USB

echo "📱 iPhone Log Access Script"
echo "=========================="

# Check if iPhone is connected
echo "🔍 Checking for connected iPhone..."
if system_profiler SPUSBDataType | grep -i "iphone\|ios" > /dev/null; then
    echo "✅ iPhone detected!"
else
    echo "❌ iPhone not detected. Please connect your iPhone via USB and trust this computer."
    echo "   Then run this script again."
    exit 1
fi

# Create output directory
OUTPUT_DIR="/Users/m/Code/iphone_logs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
echo "📁 Output directory: $OUTPUT_DIR"

# 1. Access iOS device logs via Console
echo "📋 Accessing iOS device logs..."
{
    echo "=== iOS DEVICE LOGS ==="
    echo "Date: $(date)"
    echo ""
    
    # Get device info
    echo "=== DEVICE INFO ==="
    system_profiler SPUSBDataType | grep -A 20 -i "iphone\|ios"
    echo ""
    
    # Access notification logs
    echo "=== NOTIFICATION LOGS (Last 30 days) ==="
    log show --device --predicate 'category == "notification" OR category == "NotificationCenter" OR category == "UNNotificationService"' --last 30d --style compact
    echo ""
    
    # Access message logs
    echo "=== MESSAGE LOGS (Last 30 days) ==="
    log show --device --predicate 'eventMessage CONTAINS "message" OR eventMessage CONTAINS "Message" OR category CONTAINS "Message"' --last 30d --style compact
    echo ""
    
    # Access app-specific logs
    echo "=== APP-SPECIFIC LOGS (Last 30 days) ==="
    for app in "Messages" "WhatsApp" "Telegram" "Signal" "Slack"; do
        echo "--- $app ---"
        log show --device --predicate "process == '$app'" --last 30d --style compact | head -10
        echo ""
    done
    
    # Access system logs
    echo "=== SYSTEM LOGS (Last 30 days) ==="
    log show --device --predicate 'category == "system" OR category == "System"' --last 30d --style compact | head -20
    echo ""
    
} > "$OUTPUT_DIR/iphone_logs.txt"

# 2. Access Console.app logs
echo "🖥️ Accessing Console.app logs..."
{
    echo "=== CONSOLE LOGS ==="
    echo "Date: $(date)"
    echo ""
    
    # Notification center logs
    echo "=== NOTIFICATION CENTER ==="
    log show --predicate 'subsystem == "com.apple.notificationcenter"' --last 30d --style compact
    echo ""
    
    # Focus mode logs
    echo "=== FOCUS MODES ==="
    log show --predicate 'subsystem == "com.apple.focus"' --last 30d --style compact
    echo ""
    
    # Do Not Disturb logs
    echo "=== DO NOT DISTURB ==="
    log show --predicate 'subsystem == "com.apple.DoNotDisturb"' --last 30d --style compact
    echo ""
    
} > "$OUTPUT_DIR/console_logs.txt"

# 3. Access device crash logs
echo "💥 Accessing device crash logs..."
if [ -d "$HOME/Library/Logs/CrashReporter/MobileDevice" ]; then
    {
        echo "=== DEVICE CRASH LOGS ==="
        echo "Date: $(date)"
        echo ""
        ls -la "$HOME/Library/Logs/CrashReporter/MobileDevice"
        echo ""
        echo "=== RECENT CRASH LOGS ==="
        find "$HOME/Library/Logs/CrashReporter/MobileDevice" -name "*.crash" -mtime -30 -exec ls -la {} \;
    } > "$OUTPUT_DIR/crash_logs.txt"
else
    echo "No crash logs found" > "$OUTPUT_DIR/crash_logs.txt"
fi

# 4. Access device logs via idevicesyslog (if available)
echo "📱 Checking for idevicesyslog..."
if command -v idevicesyslog >/dev/null 2>&1; then
    echo "✅ idevicesyslog found, capturing live logs..."
    timeout 30 idevicesyslog > "$OUTPUT_DIR/live_logs.txt" 2>&1 &
    LIVE_PID=$!
    sleep 30
    kill $LIVE_PID 2>/dev/null
    echo "✅ Live logs captured"
else
    echo "⚠️ idevicesyslog not found. Install with: brew install libimobiledevice"
    echo "idevicesyslog not available" > "$OUTPUT_DIR/live_logs.txt"
fi

# 5. Generate summary
echo "📄 Generating summary..."
{
    echo "=== IPHONE LOG ACCESS SUMMARY ==="
    echo "Date: $(date)"
    echo "Output Directory: $OUTPUT_DIR"
    echo ""
    echo "Files Generated:"
    ls -la "$OUTPUT_DIR"
    echo ""
    echo "=== LOG STATISTICS ==="
    echo "iPhone Logs: $(wc -l < "$OUTPUT_DIR/iphone_logs.txt") lines"
    echo "Console Logs: $(wc -l < "$OUTPUT_DIR/console_logs.txt") lines"
    echo "Crash Logs: $(wc -l < "$OUTPUT_DIR/crash_logs.txt") lines"
    echo "Live Logs: $(wc -l < "$OUTPUT_DIR/live_logs.txt") lines"
    echo ""
    echo "=== QUICK ACCESS ==="
    echo "View iPhone logs: cat $OUTPUT_DIR/iphone_logs.txt"
    echo "View Console logs: cat $OUTPUT_DIR/console_logs.txt"
    echo "View Crash logs: cat $OUTPUT_DIR/crash_logs.txt"
    echo ""
    echo "=== SEARCH FOR NOTIFICATIONS ==="
    echo "grep -i 'notification' $OUTPUT_DIR/iphone_logs.txt"
    echo "grep -i 'message' $OUTPUT_DIR/iphone_logs.txt"
    echo "grep -i 'focus' $OUTPUT_DIR/console_logs.txt"
    
} > "$OUTPUT_DIR/SUMMARY.txt"

echo ""
echo "✅ iPhone log access completed!"
echo "📁 Results saved to: $OUTPUT_DIR"
echo "📄 Summary: $OUTPUT_DIR/SUMMARY.txt"
echo ""
echo "🔍 Quick search commands:"
echo "  grep -i 'notification' $OUTPUT_DIR/iphone_logs.txt"
echo "  grep -i 'message' $OUTPUT_DIR/iphone_logs.txt"
echo "  grep -i 'focus' $OUTPUT_DIR/console_logs.txt"