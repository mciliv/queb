# Queb - Frequently Asked Questions

## General Questions

### What is Queb?

Queb is an AI-powered molecular analysis app that identifies and visualizes the chemical composition of anything you can describe, photograph, or point your camera at. It transforms everyday objects into interactive 3D molecular structures.

### How does Queb work?

1. **Input**: You provide text, image, photo, or URL
2. **AI Analysis**: OpenAI GPT-4o identifies the object and its chemical components
3. **Chemical Lookup**: Names are resolved to molecular structures via PubChem
4. **3D Generation**: RDKit creates 3D coordinates for each molecule
5. **Visualization**: 3Dmol.js renders interactive molecular structures

### Is Queb free to use?

Yes! Queb is completely free to use. All features are available without any payment or subscription.

### Do I need to create an account?

No account is required for basic usage. Accounts are optional for:
- Saving analysis history
- Accessing premium features
- API usage tracking

## Technical Questions

### What is SMILES notation?

SMILES (Simplified Molecular Input Line Entry System) is a text-based way to represent molecular structures. For example:
- Water: `O`
- Ethanol: `CCO`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`

Queb converts these to 3D structures automatically.

### What file format does Queb use for 3D structures?

Queb generates SDF (Structure Data Format) files, a standard format in chemistry that contains:
- 3D atomic coordinates
- Bond information
- Molecular properties

These files are compatible with most molecular visualization software.

### Can I download the molecular structures?

Yes! You can:
1. Right-click on the 3D viewer → "Save as Image" (PNG)
2. Access the SDF files directly via the API
3. Export SMILES notation for use in other tools

### What browsers are supported?

Queb works best on:
- Chrome/Edge (recommended)
- Firefox
- Safari
- Mobile browsers with WebGL support

Requirements:
- WebGL enabled for 3D visualization
- JavaScript enabled
- Camera permissions for camera mode

## Usage Questions

### Why does it say "No molecules found"?

Common reasons:
1. **Too vague**: Try being more specific (e.g., "coffee" instead of "drink")
2. **Non-physical object**: Queb analyzes physical substances with chemical composition
3. **Typo**: Check spelling of chemical/object names
4. **Network issue**: Ensure internet connection is stable

### How accurate is the chemical analysis?

Accuracy depends on:
- **Known compounds**: Very accurate (uses PubChem database)
- **Common objects**: Generally accurate (trained on extensive data)
- **Complex mixtures**: Shows main components, may miss trace elements
- **Novel items**: Best effort based on AI understanding

### Can Queb analyze any image?

Queb works best with:
- Clear, well-lit photos
- Distinct objects
- Common materials and substances

Limitations:
- Abstract concepts have no molecular structure
- Very dark or blurry images may fail
- Microscopic structures need specialized tools

### Why is the 3D viewer not loading?

Troubleshooting steps:
1. **Check WebGL**: Visit [get.webgl.org](https://get.webgl.org/)
2. **Update browser**: Use latest version
3. **Disable extensions**: Ad blockers may interfere
4. **Clear cache**: Force refresh with Ctrl+F5 (Cmd+R on Mac)
5. **Try different browser**: Chrome usually works best

## Privacy & Security

### Is my data private?

- **No login required**: Anonymous usage
- **No permanent storage**: Images processed and discarded
- **Optional accounts**: You control what's saved
- **Encrypted connections**: HTTPS in production

### What happens to uploaded images?

1. Image sent to server for analysis
2. AI extracts object information
3. Image data is discarded after processing
4. Only text results are returned

### Are the molecules safe to view?

Yes! Queb only displays digital 3D models. You're not handling actual chemicals. The app is educational and completely safe to use.

## Advanced Questions

### Can I use Queb's API?

Yes! See the [API Documentation](docs/API.md) for:
- REST endpoints
- Authentication
- Rate limits
- Code examples

### How can I contribute?

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Development setup
- Code guidelines
- Bug reports
- Feature requests

### Can Queb analyze proteins or DNA?

Currently, Queb focuses on small molecules. Large biomolecules like proteins and DNA are:
- Too complex for real-time rendering
- Better served by specialized tools
- On our roadmap for future versions

### Is there an offline version?

Not currently. Queb requires:
- OpenAI API for analysis
- PubChem API for chemical data
- Server processing for structure generation

### Can I integrate Queb into my app?

Future plans include:
- Embeddable widget
- SDK for developers
- White-label options

Contact us for partnership opportunities.

## Troubleshooting

### "Server error" messages

Try:
1. Refresh the page
2. Wait a moment and retry (may be rate limited)
3. Check your internet connection
4. Try a simpler query first

### Camera not working

Ensure:
1. Camera permissions granted
2. Using HTTPS (required for camera)
3. No other apps using camera
4. Browser supports getUserMedia API

### Slow performance

Optimize by:
1. Closing other browser tabs
2. Using a modern browser
3. Limiting number of molecules displayed
4. Disabling browser extensions

## Still have questions?

- **GitHub Issues**: [Report bugs or request features](https://github.com/mciliv/queb/issues)
- **Discussions**: [Ask the community](https://github.com/mciliv/queb/discussions)
- **Documentation**: Browse our [detailed docs](docs/)

---

*Last updated: October 2024*
