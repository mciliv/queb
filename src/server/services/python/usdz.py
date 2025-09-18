import os
from openbabel import openbabel
from pxr import Usd, UsdGeom, UsdShade, Sdf

def smiles_to_usdz(smiles: str, output_usdz: str):
    # Step 1: Convert SMILES to 3D coordinates using Open Babel
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", "mol")
    
    molecule = openbabel.OBMol()
    ob_conversion.ReadString(molecule, smiles)
    
    # Generate 3D coordinates
    builder = openbabel.OBBuilder()
    builder.Build(molecule)
    
    # Save as MOL file
    mol_filename = "temp.mol"
    ob_conversion.WriteFile(molecule, mol_filename)
    
    # Step 2: Export MOL to OBJ (or another 3D format)
    obj_filename = "temp.obj"
    os.system(f"obabel {mol_filename} -O {obj_filename}")
    
    # Step 3: Convert OBJ to USDZ
    stage = Usd.Stage.CreateNew(output_usdz)
    root = UsdGeom.Xform.Define(stage, "/Root")
    
    # Create a mesh for the OBJ file
    mesh = UsdGeom.Mesh.Define(stage, "/Root/Molecule")
    mesh.GetPrim().GetReferences().AddReference(obj_filename)
    
    # Optionally add a simple material
    material = UsdShade.Material.Define(stage, "/Root/Material")
    shader = UsdShade.Shader.Define(stage, "/Root/Material/Shader")
    shader.CreateIdAttr("UsdPreviewSurface")
    material.CreateSurfaceOutput().ConnectToSource(shader, "out")
    mesh.GetMaterialBind().ConnectToSource(material)
    
    # Save the USDZ file
    stage.GetRootLayer().Save()
    
    # Cleanup temporary files
    os.remove(mol_filename)
    os.remove(obj_filename)
