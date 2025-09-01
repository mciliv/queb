# Estimate & visualize what molecules are contained in a given section(s) of space

describe an object/take a photo/Upload an image → Get interactive 3D molecules

**Example**: Photo of coffee → Identifies caffeine, water, etc. → Displays 3D molecular models

UI structure:
# Text input
# camera|image|link input buttons in which 0/1 can be selected
## Camera
    Click area in the image in which the click position is the center of a red square, whose crop is used as input for identification
# Molecule x object grid
upon entered input, 
Options:
1/Single colmn
2/a column is added to the right if a column already exists
A column contains a header, specified by the object specification or inference (Load immediately upon input submittance)
each molecule have name title which links to wiki 
3dmol.js gridviewer seems best

OPEN_API_KEY is set in the .env file, but is intended to be hidden from agents via .cursorignore