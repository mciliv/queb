import React, { useEffect, useRef, useState } from 'react';
import { Terminal } from '@xterm/xterm';
import '@xterm/xterm/css/xterm.css';

const TerminalSection = ({ onCommand, isProcessing }) => {
  const terminalRef = useRef(null);
  const xtermRef = useRef(null);
  const commandRef = useRef('');

  useEffect(() => {
    if (!terminalRef.current) return;

    // Initialize xterm
    const term = new Terminal({
      cursorBlink: true,
      theme: {
        background: '#1e1e1e',
        foreground: '#f0f0f0',
        cursor: '#00ff00',
        selectionBackground: 'rgba(0, 255, 0, 0.3)'
      },
      fontFamily: 'Menlo, Monaco, "Courier New", monospace',
      fontSize: 14,
      rows: 24,
      cols: 80
    });

    term.open(terminalRef.current);
    term.write('Welcome to 1 CLI v1.0.0\r\n');
    term.write('Type a chemical name, object, or command.\r\n');
    term.write('\x1b[32m$\x1b[0m '); // Green prompt

    xtermRef.current = term;

    // Handle input
    term.onData(data => {
      const code = data.charCodeAt(0);

      // Enter key
      if (code === 13) {
        term.write('\r\n');
        const command = commandRef.current.trim();
        
        if (command) {
          handleCommand(command, term);
        } else {
          term.write('\x1b[32m$\x1b[0m ');
        }
        
        commandRef.current = '';
      } 
      // Backspace
      else if (code === 127) {
        if (commandRef.current.length > 0) {
          commandRef.current = commandRef.current.slice(0, -1);
          term.write('\b \b');
        }
      } 
      // Control characters (ignore mostly)
      else if (code < 32) {
        return;
      }
      // Normal characters
      else {
        commandRef.current += data;
        term.write(data);
      }
    });

    return () => {
      term.dispose();
    };
  }, []);

  // Handle external processing state
  useEffect(() => {
    if (!xtermRef.current) return;
    const term = xtermRef.current;

    if (isProcessing) {
      term.write('\x1b[33mProcessing...\x1b[0m\r\n');
    } else {
      // Restore prompt after processing
      term.write('\x1b[32m$\x1b[0m '); 
    }
  }, [isProcessing]);

  const handleCommand = (cmd, term) => {
    switch (cmd.toLowerCase()) {
      case 'clear':
        term.clear();
        term.write('\x1b[32m$\x1b[0m ');
        break;
      case 'help':
        term.write('Commands:\r\n');
        term.write('  [name]   - Analyze object or chemical (e.g. "coffee")\r\n');
        term.write('  clear    - Clear terminal\r\n');
        term.write('  help     - Show this help\r\n');
        term.write('\x1b[32m$\x1b[0m ');
        break;
      default:
        // Pass to parent component (App.jsx)
        onCommand(cmd);
        // Prompt will be restored after processing
        break;
    }
  };

  return (
    <div className="terminal-container" style={{ 
      padding: '10px', 
      background: '#000', 
      borderRadius: '8px',
      border: '1px solid #333',
      height: '400px',
      overflow: 'hidden'
    }}>
      <div ref={terminalRef} style={{ height: '100%', width: '100%' }} />
    </div>
  );
};

export default TerminalSection;
