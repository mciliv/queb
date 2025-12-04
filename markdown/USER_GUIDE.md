# User Guide

## Input Modes

| Mode | Icon | Use Case |
|------|------|----------|
| Text | 📝 | Type object name |
| Camera | 📷 | Live camera analysis |
| Photo | 📸 | Upload image |
| Link | 🔗 | Analyze web page |

## 3D Controls

| Action | Mouse | Touch |
|--------|-------|-------|
| Rotate | Left drag | One finger |
| Zoom | Scroll | Pinch |
| Pan | Right drag | Two fingers |
| Reset | Double click | — |

## Atom Colors (CPK)

| Color | Element |
|-------|---------|
| White | Hydrogen |
| Black/Gray | Carbon |
| Red | Oxygen |
| Blue | Nitrogen |
| Yellow | Sulfur |
| Green | Chlorine |

## API Endpoints

```
POST /api/structuralize
{ "text": "caffeine", "lookupMode": "database" }

POST /api/generate-sdfs
{ "smiles": ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"] }
```

## Tips

- Be specific: "ibuprofen tablet" not "medicine"
- Use common names: "aspirin" or "acetylsalicylic acid"
- Describe context: "stomach acid" not "acid"
