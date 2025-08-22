# Estimate & visualize what molecules are contained in a given section(s) of space

Upload an image, take a photo, or describe an object → Get molecules → See interactive 3D molecular structures

**Example**: Photo of coffee → Identifies caffeine, water, etc. → Displays 3D molecular models

UI structure:
# Text input (⌘k)
# camera|image|link input buttons in which 0/1 can be selected
## Camera 
### Mobile
    Rectangle reticle in center of screen, click screen to capture
### Non-mobile
    Click area in the image in which the click position is the center of a red square, whose crop is used as input for identification

# Molecule x object grid
upon entered input, 
Options:
1/Single colmn
2/a column is added to the right if a column already exists
A column contains a header, specified by the object specification or inference (Load immediately upon input submittance)
each molecule have name title which links to wiki 
3dmoljs gridviewer seems best

App to run on `https://localhost:3000`

## Project Structure

### Focused Codebase
- **Separated Research Code**: Molecular docking research archived separately
- **Single Purpose**: Web app focused solely on molecular space analysis
- **Clear Dependencies**: Only essential chemistry conversion tools included

### Modular Design
- **Service Layer**: Separate business logic from API routes
- **Test Organization**: Tests grouped by purpose and scope

### Performance
- **Caching**: SDF files cached to avoid regeneration
- **Optimized Rendering**: Efficient 3D molecular display

## License

MIT License - See LICENSE file for details.