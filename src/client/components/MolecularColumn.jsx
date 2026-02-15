import React from 'react';
import MoleculeViewer from './MoleculeViewer';

const MolecularColumn = ({ column, onRemove, showRemove = true }) => {
  // Render column header and viewers
  
  return (
    <div className="column">
      {/* GUARANTEED VISIBLE HEADER */}
      <div className="column-header">
        <div className="column-meta">
          <div className="column-title">
            {column.query}
          </div>
        </div>
        {showRemove && (
          <button 
            onClick={onRemove}
            className="btn-ghost"
          >
            ×
          </button>
        )}
      </div>

      {/* Failure indicator */}
      {column.failed && (
        <div className="alert-warning">
          💡 AI failed to analyze - please report this
        </div>
      )}

      {/* No molecules found (not a failure) */}
      {!column.loading && !column.failed && column.viewers.length === 0 && (
        <div className="alert-info">
          No specific molecules found. This is normal for concepts like "{column.query}".
        </div>
      )}

      {false && column.loading && column.viewers.length === 0 && (
        <div className="analyzing">
          Analyzing molecules...
        </div>
      )}

      {column.viewers.map((mol, idx) => (
        <MoleculeViewer key={idx} molecularData={mol} />
      ))}
    </div>
  );
};

export default MolecularColumn;