#!/usr/bin/env python3
"""
Clipboard Testing Utility

This script provides functionality to test clipboard operations:
- Read from clipboard
- Write to clipboard
- Test various data types
"""

import sys
import subprocess
import json
from typing import Any, Optional


class ClipboardTester:
    """A utility class for testing clipboard operations on macOS."""
    
    def __init__(self):
        self.platform = sys.platform
        if self.platform != 'darwin':
            print("Warning: This utility is optimized for macOS (darwin)")
    
    def read_clipboard(self) -> str:
        """Read text from the clipboard."""
        try:
            result = subprocess.run(['pbpaste'], capture_output=True, text=True, check=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error reading clipboard: {e}")
            return ""
        except FileNotFoundError:
            print("Error: pbpaste command not found. Make sure you're on macOS.")
            return ""
    
    def write_clipboard(self, text: str) -> bool:
        """Write text to the clipboard."""
        try:
            subprocess.run(['pbcopy'], input=text, text=True, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error writing to clipboard: {e}")
            return False
        except FileNotFoundError:
            print("Error: pbcopy command not found. Make sure you're on macOS.")
            return False
    
    def test_clipboard_read(self) -> dict:
        """Test reading from clipboard and return analysis."""
        print("Testing clipboard read...")
        content = self.read_clipboard()
        
        result = {
            'success': bool(content),
            'content': content,
            'length': len(content),
            'type': type(content).__name__,
            'is_empty': len(content) == 0,
            'lines': content.count('\n') + 1 if content else 0,
            'words': len(content.split()) if content else 0
        }
        
        print(f"✓ Clipboard read successful: {result['success']}")
        print(f"  Content length: {result['length']} characters")
        print(f"  Lines: {result['lines']}")
        print(f"  Words: {result['words']}")
        
        return result
    
    def test_clipboard_write(self, test_data: str = "Hello from clipboard test!") -> dict:
        """Test writing to clipboard."""
        print(f"Testing clipboard write with: '{test_data}'")
        
        success = self.write_clipboard(test_data)
        
        result = {
            'success': success,
            'test_data': test_data,
            'data_length': len(test_data)
        }
        
        if success:
            print("✓ Clipboard write successful")
            # Verify by reading back
            read_back = self.read_clipboard()
            result['verified'] = read_back == test_data
            print(f"✓ Verification: {'PASSED' if result['verified'] else 'FAILED'}")
        else:
            print("✗ Clipboard write failed")
        
        return result
    
    def test_various_data_types(self):
        """Test clipboard with various data types."""
        test_cases = [
            ("Simple text", "Hello World"),
            ("Multiline text", "Line 1\nLine 2\nLine 3"),
            ("JSON data", json.dumps({"name": "test", "value": 123})),
            ("Empty string", ""),
            ("Special characters", "!@#$%^&*()_+-=[]{}|;':\",./<>?"),
            ("Unicode text", "Hello 世界 🌍"),
            ("Large text", "A" * 1000)
        ]
        
        print("\nTesting various data types...")
        results = []
        
        for name, data in test_cases:
            print(f"\nTesting: {name}")
            write_result = self.test_clipboard_write(data)
            read_result = self.test_clipboard_read()
            
            results.append({
                'name': name,
                'write_success': write_result['success'],
                'read_success': read_result['success'],
                'data_preserved': read_result['content'] == data,
                'data_length': len(data)
            })
        
        return results
    
    def interactive_mode(self):
        """Run interactive clipboard testing mode."""
        print("=== Clipboard Testing Utility ===")
        print("Commands:")
        print("  read    - Read from clipboard")
        print("  write   - Write to clipboard")
        print("  test    - Run all tests")
        print("  quit    - Exit")
        print()
        
        while True:
            try:
                command = input("Enter command: ").strip().lower()
                
                if command == 'quit':
                    break
                elif command == 'read':
                    self.test_clipboard_read()
                elif command == 'write':
                    text = input("Enter text to write: ")
                    self.test_clipboard_write(text)
                elif command == 'test':
                    self.test_various_data_types()
                else:
                    print("Unknown command. Try: read, write, test, quit")
                
                print()
                
            except KeyboardInterrupt:
                print("\nExiting...")
                break


def main():
    """Main function to run clipboard tests."""
    tester = ClipboardTester()
    
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == 'read':
            result = tester.test_clipboard_read()
            print(f"\nCurrent clipboard content:\n{result['content']}")
        elif command == 'write' and len(sys.argv) > 2:
            text = ' '.join(sys.argv[2:])
            tester.test_clipboard_write(text)
        elif command == 'test':
            tester.test_various_data_types()
        elif command == 'interactive':
            tester.interactive_mode()
        else:
            print("Usage: python clipboard_test.py [read|write <text>|test|interactive]")
    else:
        # Default: run basic tests
        print("Running basic clipboard tests...")
        tester.test_clipboard_read()
        print()
        tester.test_clipboard_write("Test data from clipboard utility")
        print()
        print("Run 'python clipboard_test.py interactive' for interactive mode")


if __name__ == "__main__":
    main()