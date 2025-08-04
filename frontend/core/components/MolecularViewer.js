// MolecularViewer.js - 3D molecular visualization and SDF processing
export class MolecularViewer {
  constructor() {
    this.viewers = [];
    this.snapshots = null;
  }

  async initialize() {
    this.snapshots = document.querySelector(".snapshots-container");
    this.setupEventListeners();
    console.log('✅ Molecular viewer initialized');
  }

  setupEventListeners() {
    // Listen for analysis results
    document.addEventListener('analysisResult', async (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      await this.handleAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });
  }

  async handleAnalysisResult(output, icon, objectName, useQuotes, croppedImageData) {
    try {
      const chemicals = this.extractChemicals(output);
      const sdfFiles = await this.generateSDFs(chemicals);
      
      await this.createObjectColumn(
        objectName, 
        sdfFiles, 
        chemicals.smiles, 
        null, 
        chemicals.summary, 
        chemicals.skippedChemicals, 
        output.description, 
        chemicals.chemicals, 
        croppedImageData
      );
    } catch (error) {
      console.error('❌ Error processing analysis result:', error);
      document.dispatchEvent(new CustomEvent('appError', { 
        detail: { message: 'Failed to process analysis result', type: 'error' } 
      }));
    }
  }

  extractChemicals(output) {
    const chemicals = [];
    const smiles = [];
    const skippedChemicals = [];
    let summary = null;
    let description = null;

    if (output.chemicals) {
      output.chemicals.forEach(chem => {
        if (chem.smiles) {
          chemicals.push(chem);
          smiles.push(chem.smiles);
        } else {
          skippedChemicals.push(chem);
        }
      });
    }

    if (output.summary) {
      summary = output.summary;
    }

    if (output.description) {
      description = output.description;
    }

    return { chemicals, smiles, skippedChemicals, summary, description };
  }

  async generateSDFs(chemicals) {
    const sdfFiles = [];
    
    for (const smiles of chemicals.smiles) {
      try {
        const response = await fetch('/api/generate-sdf', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles })
        });

        if (!response.ok) throw new Error(`SDF generation failed for ${smiles}`);
        
        const sdfData = await response.text();
        sdfFiles.push({ smiles, sdf: sdfData });
      } catch (error) {
        console.error(`❌ Failed to generate SDF for ${smiles}:`, error);
      }
    }

    return sdfFiles;
  }

  async createObjectColumn(objectName, sdfFiles, smiles = [], errorMessage = null, 
                          summary = null, skippedChemicals = [], description = null, 
                          chemicals = null, croppedImageData = null) {
    const column = document.createElement('div');
    column.className = 'object-column';
    column.style.cssText = `
      display: flex;
      flex-direction: column;
      gap: 10px;
      padding: 15px;
      background: rgba(255, 255, 255, 0.05);
      border-radius: 8px;
      margin-bottom: 15px;
    `;

    // Header with close button
    const header = document.createElement('div');
    header.style.cssText = `
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 10px;
    `;

    const title = document.createElement('h3');
    title.textContent = objectName;
    title.style.cssText = `
      margin: 0;
      font-size: 16px;
      font-weight: 600;
      color: #00d4ff;
    `;

    const closeBtn = document.createElement('button');
    closeBtn.innerHTML = '✕';
    closeBtn.style.cssText = `
      background: none;
      border: none;
      color: #666;
      cursor: pointer;
      font-size: 18px;
      padding: 5px;
      border-radius: 4px;
      transition: color 0.2s;
    `;
    closeBtn.onmouseover = () => closeBtn.style.color = '#ff4444';
    closeBtn.onmouseout = () => closeBtn.style.color = '#666';
    closeBtn.onclick = () => column.remove();

    header.appendChild(title);
    header.appendChild(closeBtn);
    column.appendChild(header);

    // Error message if any
    if (errorMessage) {
      const errorDiv = document.createElement('div');
      errorDiv.textContent = errorMessage;
      errorDiv.style.cssText = `
        color: #ff4444;
        padding: 10px;
        background: rgba(255, 68, 68, 0.1);
        border-radius: 4px;
        margin-bottom: 10px;
      `;
      column.appendChild(errorDiv);
    }

    // Description if available
    if (description) {
      const descDiv = document.createElement('div');
      descDiv.textContent = description;
      descDiv.style.cssText = `
        color: #ccc;
        font-size: 14px;
        margin-bottom: 10px;
        line-height: 1.4;
      `;
      column.appendChild(descDiv);
    }

    // Summary if available
    if (summary) {
      const summaryDiv = this.createSummaryElement(summary);
      column.appendChild(summaryDiv);
    }

    // Molecular viewers
    if (sdfFiles.length > 0) {
      const viewerContainer = document.createElement('div');
      viewerContainer.style.cssText = `
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 10px;
        margin-bottom: 10px;
      `;

      await this.createGridViewer(viewerContainer, sdfFiles, chemicals, smiles);
      column.appendChild(viewerContainer);
    }

    // Skipped chemicals if any
    if (skippedChemicals.length > 0) {
      const skippedDiv = this.createSkippedChemicalsElement(skippedChemicals);
      column.appendChild(skippedDiv);
    }

    this.snapshots.appendChild(column);
    this.updateScrollHandles();
  }

  createSummaryElement(summary) {
    const summaryDiv = document.createElement('div');
    summaryDiv.className = 'chemical-summary';
    summaryDiv.style.cssText = `
      display: flex;
      flex-wrap: wrap;
      gap: 15px;
      font-size: 14px;
      color: #ccc;
      margin-bottom: 10px;
    `;

    if (summary.total) {
      const totalDiv = document.createElement('div');
      totalDiv.innerHTML = `Total chemicals found: <span style="color: #00d4ff;">${summary.total}</span>`;
      summaryDiv.appendChild(totalDiv);
    }

    if (summary.visualizable) {
      const visualDiv = document.createElement('div');
      visualDiv.innerHTML = `3D visualizable: <span style="color: #00d4ff;">${summary.visualizable}</span>`;
      summaryDiv.appendChild(visualDiv);
    }

    if (summary.skipped && summary.skipped > 0) {
      const skippedDiv = document.createElement('div');
      skippedDiv.innerHTML = `Non-SMILES formats: <span style="color: #ffaa00;">${summary.skipped}</span>`;
      summaryDiv.appendChild(skippedDiv);
    }

    if (summary.errors && summary.errors > 0) {
      const errorsDiv = document.createElement('div');
      errorsDiv.innerHTML = `Failed: <span style="color: #ff4444;">${summary.errors}</span>`;
      summaryDiv.appendChild(errorsDiv);
    }

    return summaryDiv;
  }

  createSkippedChemicalsElement(skippedChemicals) {
    const skippedDiv = document.createElement('div');
    skippedDiv.className = 'skipped-chemicals';
    skippedDiv.style.cssText = `
      margin-top: 10px;
      padding: 10px;
      background: rgba(255, 170, 0, 0.1);
      border-radius: 4px;
    `;

    const title = document.createElement('div');
    title.textContent = 'Other chemicals found:';
    title.style.cssText = `
      font-weight: 600;
      color: #ffaa00;
      margin-bottom: 5px;
    `;
    skippedDiv.appendChild(title);

    const list = document.createElement('div');
    list.className = 'skipped-list';
    list.style.cssText = `
      font-size: 14px;
      color: #ccc;
      margin-bottom: 5px;
    `;
    
    skippedChemicals.forEach(chem => {
      const item = document.createElement('div');
      item.textContent = chem.name || chem;
      list.appendChild(item);
    });
    skippedDiv.appendChild(list);

    const note = document.createElement('div');
    note.textContent = '* These likely represent minerals/crystals that can\'t be shown in 3D molecular view';
    note.style.cssText = `
      font-size: 12px;
      color: #888;
      font-style: italic;
    `;
    skippedDiv.appendChild(note);

    return skippedDiv;
  }

  async createGridViewer(container, sdfFiles, chemicals, smiles) {
    for (let i = 0; i < sdfFiles.length; i++) {
      const sdfFile = sdfFiles[i];
      const chemical = chemicals ? chemicals[i] : null;
      const smilesStr = smiles[i];

      const viewerDiv = document.createElement('div');
      viewerDiv.style.cssText = `
        position: relative;
        background: rgba(0, 0, 0, 0.3);
        border-radius: 6px;
        overflow: hidden;
        min-height: 200px;
      `;

      const label = document.createElement('div');
      label.textContent = chemical ? chemical.name : smilesStr;
      label.style.cssText = `
        position: absolute;
        top: 10px;
        left: 10px;
        background: rgba(0, 0, 0, 0.7);
        color: #00d4ff;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 12px;
        font-weight: 500;
        z-index: 10;
        max-width: calc(100% - 20px);
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      `;
      viewerDiv.appendChild(label);

      const renderDiv = document.createElement('div');
      renderDiv.style.cssText = `
        width: 100%;
        height: 200px;
        background: transparent;
      `;
      viewerDiv.appendChild(renderDiv);

      container.appendChild(viewerDiv);

      try {
        await this.render(sdfFile.sdf, renderDiv);
      } catch (error) {
        console.error('❌ Failed to render molecule:', error);
        renderDiv.innerHTML = '<div style="color: #ff4444; padding: 20px; text-align: center;">Render failed</div>';
      }
    }
  }

  async render(sdfData, container) {
    return new Promise((resolve, reject) => {
      try {
        const viewer = $3Dmol.createViewer(container, {
          backgroundColor: 'transparent'
        });

        viewer.addModel(sdfData, 'sdf');
        viewer.setStyle({}, { sphere: { radius: 0.8 } });
        viewer.zoomTo();
        viewer.render();

        this.viewers.push(viewer);
        resolve();
      } catch (error) {
        reject(error);
      }
    });
  }

  updateScrollHandles() {
    // Update scroll handles if needed
    const scrollContainer = this.snapshots.closest('.scroll-container');
    if (scrollContainer) {
      // Trigger scroll update if scroll handles exist
      const event = new Event('scroll');
      scrollContainer.dispatchEvent(event);
    }
  }

  clearResults() {
    if (this.snapshots) {
      this.snapshots.innerHTML = '';
    }
    const gldiv = document.getElementById('gldiv');
    if (gldiv) {
      gldiv.innerHTML = '';
    }
    this.viewers = [];
  }
} 