#!/bin/bash

# iPhone Notification Center Investigation Script
# This script provides additional system-level access for investigating
# messages that only appeared in notification center this past month

echo "🚀 iPhone Notification Investigation - System Access"
echo "=================================================="

# Create output directory
OUTPUT_DIR="/Users/m/Code/iphone_investigation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
echo "📁 Output directory: $OUTPUT_DIR"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# 1. System Information
log "🔍 Gathering system information..."
{
    echo "=== SYSTEM INFORMATION ==="
    echo "Date: $(date)"
    echo "Uptime: $(uptime)"
    echo "macOS Version:"
    sw_vers
    echo ""
    echo "Hardware:"
    system_profiler SPHardwareDataType
    echo ""
    echo "Memory Usage:"
    vm_stat
    echo ""
    echo "Disk Usage:"
    df -h
} > "$OUTPUT_DIR/system_info.txt"

# 2. Check for connected iOS devices
log "📱 Checking for connected iOS devices..."
{
    echo "=== CONNECTED iOS DEVICES ==="
    system_profiler SPUSBDataType | grep -A 10 -B 2 -i "iphone\|ipad\|ios"
    echo ""
    echo "=== iOS DEVICE LOGS ==="
    if [ -d "$HOME/Library/Logs/CrashReporter/MobileDevice" ]; then
        echo "iOS device logs found:"
        ls -la "$HOME/Library/Logs/CrashReporter/MobileDevice"
    else
        echo "No iOS device logs found"
    fi
} > "$OUTPUT_DIR/ios_devices.txt"

# 3. Notification Center Logs
log "📋 Analyzing notification center logs..."
{
    echo "=== NOTIFICATION CENTER LOGS (Last 30 days) ==="
    log show --predicate 'category == "notification" OR category == "NotificationCenter" OR category == "UNNotificationService"' --last 30d --style compact
} > "$OUTPUT_DIR/notification_logs.txt"

# 4. Message-related logs
log "💬 Analyzing message-related logs..."
{
    echo "=== MESSAGE-RELATED LOGS (Last 30 days) ==="
    log show --predicate 'eventMessage CONTAINS "message" OR eventMessage CONTAINS "Message" OR category CONTAINS "Message"' --last 30d --style compact
} > "$OUTPUT_DIR/message_logs.txt"

# 5. Recent system changes
log "🔄 Checking recent system changes..."
{
    echo "=== RECENT SOFTWARE UPDATES ==="
    softwareupdate --history
    echo ""
    echo "=== RECENT APP MODIFICATIONS ==="
    find /Applications -name "*.app" -mtime -30 -exec ls -la {} \;
    echo ""
    echo "=== RECENT SYSTEM PREFERENCES CHANGES ==="
    log show --predicate 'category == "preferences" OR category == "Preferences"' --last 30d --style compact
} > "$OUTPUT_DIR/recent_changes.txt"

# 6. Focus and Do Not Disturb settings
log "🔕 Checking Focus and Do Not Disturb settings..."
{
    echo "=== FOCUS MODE LOGS ==="
    log show --predicate 'category == "Focus" OR category == "DoNotDisturb"' --last 30d --style compact
    echo ""
    echo "=== NOTIFICATION SETTINGS ==="
    defaults read com.apple.notificationcenter
} > "$OUTPUT_DIR/focus_settings.txt"

# 7. App-specific notification logs
log "📱 Analyzing app-specific notification logs..."
{
    echo "=== APP-SPECIFIC NOTIFICATION LOGS ==="
    # Common messaging apps
    for app in "Messages" "WhatsApp" "Telegram" "Signal" "Slack" "Discord"; do
        echo "--- $app ---"
        log show --predicate "process == '$app'" --last 30d --style compact | head -20
        echo ""
    done
} > "$OUTPUT_DIR/app_notifications.txt"

# 8. System performance logs
log "⚡ Checking system performance logs..."
{
    echo "=== SYSTEM PERFORMANCE LOGS ==="
    log show --predicate 'category == "performance" OR category == "Performance"' --last 30d --style compact
    echo ""
    echo "=== MEMORY PRESSURE LOGS ==="
    log show --predicate 'eventMessage CONTAINS "memory" OR eventMessage CONTAINS "Memory"' --last 30d --style compact
} > "$OUTPUT_DIR/performance_logs.txt"

# 9. Generate summary report
log "📄 Generating summary report..."
{
    echo "=== INVESTIGATION SUMMARY ==="
    echo "Investigation Date: $(date)"
    echo "Output Directory: $OUTPUT_DIR"
    echo ""
    echo "Files Generated:"
    ls -la "$OUTPUT_DIR"
    echo ""
    echo "=== KEY FINDINGS ==="
    echo "Notification Center Logs: $(wc -l < "$OUTPUT_DIR/notification_logs.txt") lines"
    echo "Message Logs: $(wc -l < "$OUTPUT_DIR/message_logs.txt") lines"
    echo "Recent Changes: $(wc -l < "$OUTPUT_DIR/recent_changes.txt") lines"
    echo "Focus Settings: $(wc -l < "$OUTPUT_DIR/focus_settings.txt") lines"
    echo "App Notifications: $(wc -l < "$OUTPUT_DIR/app_notifications.txt") lines"
    echo ""
    echo "=== NEXT STEPS ==="
    echo "1. Review notification_logs.txt for NotificationCenter entries"
    echo "2. Check message_logs.txt for message-specific patterns"
    echo "3. Examine recent_changes.txt for system updates"
    echo "4. Look at focus_settings.txt for Do Not Disturb changes"
    echo "5. Check app_notifications.txt for app-specific issues"
    echo ""
    echo "=== COMMON CAUSES OF NOTIFICATION-ONLY MESSAGES ==="
    echo "- Focus modes or Do Not Disturb settings"
    echo "- App-specific notification settings"
    echo "- Recent iOS or app updates"
    echo "- Notification grouping settings"
    echo "- Background app refresh settings"
    echo "- Notification delivery preferences"
} > "$OUTPUT_DIR/SUMMARY.txt"

# 10. Create analysis script
log "🔧 Creating analysis script..."
cat > "$OUTPUT_DIR/analyze_results.py" << 'EOF'
#!/usr/bin/env python3
"""
Analysis script for iPhone notification investigation results
"""

import os
import re
from datetime import datetime, timedelta

def analyze_logs():
    """Analyze the collected log files"""
    output_dir = os.path.dirname(os.path.abspath(__file__))
    
    print("🔍 Analyzing iPhone notification investigation results...")
    print("=" * 60)
    
    # Analyze notification center logs
    notification_file = os.path.join(output_dir, "notification_logs.txt")
    if os.path.exists(notification_file):
        with open(notification_file, 'r') as f:
            content = f.read()
            
        notification_center_count = len(re.findall(r'NotificationCenter', content))
        notification_count = len(re.findall(r'notification', content, re.IGNORECASE))
        
        print(f"📋 Notification Center entries: {notification_center_count}")
        print(f"📋 Total notification entries: {notification_count}")
    
    # Analyze message logs
    message_file = os.path.join(output_dir, "message_logs.txt")
    if os.path.exists(message_file):
        with open(message_file, 'r') as f:
            content = f.read()
            
        message_count = len(re.findall(r'message', content, re.IGNORECASE))
        print(f"💬 Message-related entries: {message_count}")
    
    # Look for patterns
    print("\n🔍 Looking for patterns...")
    
    # Check for recent changes
    changes_file = os.path.join(output_dir, "recent_changes.txt")
    if os.path.exists(changes_file):
        with open(changes_file, 'r') as f:
            content = f.read()
            
        if "softwareupdate" in content.lower():
            print("⚠️ Recent software updates detected")
        if "Applications" in content:
            print("⚠️ Recent app modifications detected")
    
    print("\n✅ Analysis complete!")
    print("📁 Check individual log files for detailed information")

if __name__ == "__main__":
    analyze_logs()
EOF

chmod +x "$OUTPUT_DIR/analyze_results.py"

log "✅ Investigation completed!"
echo ""
echo "📁 Results saved to: $OUTPUT_DIR"
echo "📄 Summary report: $OUTPUT_DIR/SUMMARY.txt"
echo "🔧 Analysis script: $OUTPUT_DIR/analyze_results.py"
echo ""
echo "💡 To analyze results, run:"
echo "   cd $OUTPUT_DIR && python3 analyze_results.py"
echo ""
echo "🔍 To view specific logs:"
echo "   cat $OUTPUT_DIR/notification_logs.txt"
echo "   cat $OUTPUT_DIR/message_logs.txt"