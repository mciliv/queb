import React, { useMemo, useState } from 'react';
import { findElement } from '../data/element-data.js';

const SHELL_CAPACITIES = [2, 8, 18, 32, 32, 18, 8];

const buildShells = (atomicNumber) => {
  let remaining = atomicNumber;
  const shells = [];
  for (const capacity of SHELL_CAPACITIES) {
    if (remaining <= 0) break;
    const count = Math.min(remaining, capacity);
    shells.push(count);
    remaining -= count;
  }
  return shells;
};

const buildElectronPositions = (electronCount, radius) => {
  // User prompt: build a atomic orbital visualizer that takes a specification of the element, eg "aluminum"
  if (electronCount === 0) return [];
  const angleStep = (Math.PI * 2) / electronCount;
  return Array.from({ length: electronCount }).map((_, index) => {
    const angle = angleStep * index;
    return {
      x: Math.cos(angle) * radius,
      y: Math.sin(angle) * radius
    };
  });
};

const AtomicOrbitalVisualizer = ({ defaultSpec = 'aluminum' }) => {
  const [input, setInput] = useState(defaultSpec);
  const [submittedSpec, setSubmittedSpec] = useState(defaultSpec);
  const [error, setError] = useState('');

  const element = useMemo(() => findElement(submittedSpec), [submittedSpec]);
  const shells = useMemo(() => (element ? buildShells(element.atomicNumber) : []), [element]);
  const electronTotal = element ? element.atomicNumber : 0;
  const shellConfig = shells.join('-');

  const handleSubmit = (event) => {
    event.preventDefault();
    const value = input.trim();
    if (!value) {
      setError('Enter an element name, symbol, or atomic number.');
      return;
    }
    const resolved = findElement(value);
    if (!resolved) {
      setError(`No element found for "${value}". Try "Al", "Aluminum", or "13".`);
      return;
    }
    setError('');
    setSubmittedSpec(value);
  };

  return (
    <div className="orbital-card">
      <div className="orbital-header">
        <div>
          <div className="orbital-title">Atomic Orbital Visualizer</div>
          <div className="orbital-subtitle">Bohr-style shell approximation</div>
        </div>
        {element && (
          <div className="orbital-chip">
            {element.name} · {element.symbol}
          </div>
        )}
      </div>

      <form className="orbital-input" onSubmit={handleSubmit}>
        <input
          type="text"
          value={input}
          onChange={(event) => {
            setInput(event.target.value);
            if (error) setError('');
          }}
          placeholder="Element name, symbol, or atomic number"
          className="input-base"
        />
        <button className="btn-icon" type="submit" aria-label="Visualize">
          →
        </button>
      </form>

      {error && <div className="error-text">{error}</div>}

      {element && (
        <div className="orbital-body">
          <div className="orbital-diagram" role="img" aria-label={`Electron shells for ${element.name}`}>
            <div className="orbital-nucleus">
              <div className="orbital-symbol">{element.symbol}</div>
              <div className="orbital-atomic-number">{element.atomicNumber}</div>
            </div>
            {shells.map((count, shellIndex) => {
              const shellRadius = 50 + shellIndex * 36;
              const positions = buildElectronPositions(count, shellRadius);
              return (
                <div
                  key={`shell-${shellIndex}`}
                  className="orbital-shell"
                  style={{
                    width: shellRadius * 2,
                    height: shellRadius * 2
                  }}
                >
                  {positions.map((pos, idx) => (
                    <div
                      key={`electron-${shellIndex}-${idx}`}
                      className="orbital-electron"
                      style={{
                        transform: `translate(${pos.x}px, ${pos.y}px)`
                      }}
                    />
                  ))}
                </div>
              );
            })}
          </div>

          <div className="orbital-metadata">
            <div>
              <div className="orbital-label">Electron count</div>
              <div className="orbital-value">{electronTotal}</div>
            </div>
            <div>
              <div className="orbital-label">Shells</div>
              <div className="orbital-value">{shellConfig}</div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default AtomicOrbitalVisualizer;
