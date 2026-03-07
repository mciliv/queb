#!/usr/bin/env python3
"""
iPhone Notification Center Investigation Script

This script helps investigate messages that only appeared in the notification center
on iPhone this past month by accessing system information and logs.
"""

import subprocess
import json
import os
import sys
from datetime import datetime, timedelta
import re

class iPhoneNotificationInvestigator:
    def __init__(self):
        self.results = {
            'investigation_date': datetime.now().isoformat(),
            'system_info': {},
            'notification_logs': [],
            'recent_changes': [],
            'findings': []
        }
    
    def get_system_info(self):
        """Gather basic system information"""
        print("🔍 Gathering system information...")
        
        try:
            # Get macOS version
            result = subprocess.run(['sw_vers'], capture_output=True, text=True)
            if result.returncode == 0:
                self.results['system_info']['macos_version'] = result.stdout.strip()
            
            # Get system uptime
            result = subprocess.run(['uptime'], capture_output=True, text=True)
            if result.returncode == 0:
                self.results['system_info']['uptime'] = result.stdout.strip()
            
            # Get available disk space
            result = subprocess.run(['df', '-h', '/'], capture_output=True, text=True)
            if result.returncode == 0:
                self.results['system_info']['disk_usage'] = result.stdout.strip()
                
        except Exception as e:
            print(f"❌ Error gathering system info: {e}")
    
    def check_ios_device_connection(self):
        """Check if iPhone is connected and accessible"""
        print("📱 Checking iPhone connection...")
        
        try:
            # Check if iPhone is connected via USB
            result = subprocess.run(['system_profiler', 'SPUSBDataType'], capture_output=True, text=True)
            if 'iPhone' in result.stdout or 'iOS' in result.stdout:
                self.results['system_info']['iphone_connected'] = True
                print("✅ iPhone detected via USB")
            else:
                self.results['system_info']['iphone_connected'] = False
                print("⚠️ iPhone not detected via USB")
            
            # Check for iOS device logs
            log_path = os.path.expanduser('~/Library/Logs/CrashReporter/MobileDevice')
            if os.path.exists(log_path):
                self.results['system_info']['ios_logs_available'] = True
                print("✅ iOS device logs found")
            else:
                self.results['system_info']['ios_logs_available'] = False
                print("⚠️ iOS device logs not found")
                
        except Exception as e:
            print(f"❌ Error checking iPhone connection: {e}")
    
    def analyze_notification_logs(self):
        """Analyze notification-related logs"""
        print("📋 Analyzing notification logs...")
        
        # Check Console.app logs for notification-related entries
        try:
            # Look for notification center logs
            result = subprocess.run([
                'log', 'show', '--predicate', 
                'category == "notification" OR category == "NotificationCenter" OR category == "UNNotificationService"',
                '--last', '30d',
                '--style', 'json'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                logs = []
                for line in result.stdout.strip().split('\n'):
                    if line:
                        try:
                            log_entry = json.loads(line)
                            logs.append(log_entry)
                        except json.JSONDecodeError:
                            continue
                
                self.results['notification_logs'] = logs
                print(f"✅ Found {len(logs)} notification-related log entries")
            else:
                print("⚠️ Could not access notification logs")
                
        except Exception as e:
            print(f"❌ Error analyzing notification logs: {e}")
    
    def check_recent_system_changes(self):
        """Check for recent system changes that might affect notifications"""
        print("🔄 Checking for recent system changes...")
        
        try:
            # Check for recent software updates
            result = subprocess.run([
                'softwareupdate', '--history'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                self.results['recent_changes']['software_updates'] = result.stdout.strip()
                print("✅ Software update history retrieved")
            
            # Check for recent app installations
            app_path = '/Applications'
            recent_apps = []
            for app in os.listdir(app_path):
                app_full_path = os.path.join(app_path, app)
                if os.path.isdir(app_full_path) and app.endswith('.app'):
                    stat = os.stat(app_full_path)
                    mod_time = datetime.fromtimestamp(stat.st_mtime)
                    if mod_time > datetime.now() - timedelta(days=30):
                        recent_apps.append({
                            'name': app,
                            'modified': mod_time.isoformat()
                        })
            
            self.results['recent_changes']['recent_apps'] = recent_apps
            print(f"✅ Found {len(recent_apps)} recently modified apps")
            
        except Exception as e:
            print(f"❌ Error checking recent changes: {e}")
    
    def analyze_notification_patterns(self):
        """Analyze patterns in notification data"""
        print("🔍 Analyzing notification patterns...")
        
        findings = []
        
        # Look for notification center specific entries
        notification_center_entries = [
            log for log in self.results['notification_logs']
            if 'NotificationCenter' in str(log.get('category', ''))
        ]
        
        if notification_center_entries:
            findings.append({
                'type': 'notification_center_activity',
                'count': len(notification_center_entries),
                'description': f'Found {len(notification_center_entries)} NotificationCenter entries in the last 30 days'
            })
        
        # Look for message-related notifications
        message_notifications = [
            log for log in self.results['notification_logs']
            if 'message' in str(log.get('message', '')).lower() or 'Message' in str(log.get('category', ''))
        ]
        
        if message_notifications:
            findings.append({
                'type': 'message_notifications',
                'count': len(message_notifications),
                'description': f'Found {len(message_notifications)} message-related notifications'
            })
        
        # Look for recent changes that might affect notifications
        if isinstance(self.results['recent_changes'], dict) and self.results['recent_changes'].get('recent_apps'):
            findings.append({
                'type': 'recent_app_changes',
                'count': len(self.results['recent_changes']['recent_apps']),
                'description': f'Found {len(self.results["recent_changes"]["recent_apps"])} recently modified apps that might affect notifications'
            })
        
        self.results['findings'] = findings
        
        for finding in findings:
            print(f"📊 {finding['description']}")
    
    def generate_report(self):
        """Generate investigation report"""
        print("\n📄 Generating investigation report...")
        
        report_path = f"/Users/m/Code/iphone_notification_investigation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        
        with open(report_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        print(f"✅ Report saved to: {report_path}")
        
        # Generate summary
        print("\n" + "="*60)
        print("📋 INVESTIGATION SUMMARY")
        print("="*60)
        
        print(f"Investigation Date: {self.results['investigation_date']}")
        print(f"iPhone Connected: {self.results['system_info'].get('iphone_connected', 'Unknown')}")
        print(f"iOS Logs Available: {self.results['system_info'].get('ios_logs_available', 'Unknown')}")
        print(f"Notification Log Entries: {len(self.results['notification_logs'])}")
        print(f"Recent App Changes: {len(self.results['recent_changes'].get('recent_apps', []))}")
        
        print("\n🔍 Key Findings:")
        for finding in self.results['findings']:
            print(f"  • {finding['description']}")
        
        print("\n💡 Next Steps:")
        print("  1. Check the detailed JSON report for specific log entries")
        print("  2. Look for patterns in notification timing and content")
        print("  3. Check iPhone Settings > Notifications for app-specific settings")
        print("  4. Consider checking Focus modes and Do Not Disturb settings")
        print("  5. Review recent app updates that might affect notification behavior")
    
    def run_investigation(self):
        """Run the complete investigation"""
        print("🚀 Starting iPhone Notification Investigation")
        print("="*50)
        
        self.get_system_info()
        self.check_ios_device_connection()
        self.analyze_notification_logs()
        self.check_recent_system_changes()
        self.analyze_notification_patterns()
        self.generate_report()
        
        print("\n✅ Investigation completed!")

if __name__ == "__main__":
    investigator = iPhoneNotificationInvestigator()
    investigator.run_investigation()