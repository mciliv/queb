import React, { useRef, useMemo } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';

const CPK_COLORS = {
  H: '#FFFFFF', C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050',
  Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6', Si: '#F0C8A0', P: '#FF8000',
  S: '#FFFF30', Cl: '#1FF01F', K: '#8F40D4', Ca: '#3DFF00', Fe: '#E06633',
  Cu: '#C88033', Zn: '#7D80B0', Br: '#A62929', I: '#940094', Au: '#FFD123',
  Ag: '#C0C0C0', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5'
};

const VDW_SCALE = 0.3; // Scale down van der Waals radii for visualization

const VDW_RADII = {
  H: 1.20, C: 1.70, N: 1.55, O: 1.52, F: 1.47, Na: 2.27, Cl: 1.75,
  K: 2.75, Ca: 2.31, Fe: 2.04, Cu: 1.40, S: 1.80, P: 1.80, Si: 2.10,
  Au: 1.66, Ag: 1.72, Li: 1.82, Al: 1.84, Br: 1.85, I: 1.98
};

function AtomInstances({ atoms }) {
  const meshRef = useRef();

  // Group atoms by element for instanced rendering
  const groups = useMemo(() => {
    const map = new Map();
    for (const atom of atoms) {
      const el = atom.element || 'C';
      if (!map.has(el)) map.set(el, []);
      map.get(el).push(atom);
    }
    return Array.from(map.entries());
  }, [atoms]);

  return (
    <>
      {groups.map(([element, elementAtoms]) => (
        <AtomGroup key={element} element={element} atoms={elementAtoms} />
      ))}
    </>
  );
}

function AtomGroup({ element, atoms }) {
  const meshRef = useRef();
  const color = CPK_COLORS[element] || '#FF69B4';
  const radius = (VDW_RADII[element] || 1.5) * VDW_SCALE;

  const geometry = useMemo(() => new THREE.SphereGeometry(1, 16, 12), []);
  const material = useMemo(() => new THREE.MeshStandardMaterial({
    color: new THREE.Color(color),
    roughness: 0.4,
    metalness: 0.1
  }), [color]);

  const matrices = useMemo(() => {
    const temp = new THREE.Matrix4();
    const arr = new Float32Array(atoms.length * 16);
    atoms.forEach((atom, i) => {
      temp.makeTranslation(atom.x, atom.y, atom.z);
      temp.scale(new THREE.Vector3(radius, radius, radius));
      temp.toArray(arr, i * 16);
    });
    return arr;
  }, [atoms, radius]);

  React.useEffect(() => {
    if (!meshRef.current) return;
    const mesh = meshRef.current;
    const temp = new THREE.Matrix4();
    atoms.forEach((atom, i) => {
      temp.makeTranslation(atom.x, atom.y, atom.z);
      temp.scale(new THREE.Vector3(radius, radius, radius));
      mesh.setMatrixAt(i, temp);
    });
    mesh.instanceMatrix.needsUpdate = true;
  }, [atoms, radius]);

  return (
    <instancedMesh ref={meshRef} args={[geometry, material, atoms.length]} />
  );
}

function SceneContent({ sceneData }) {
  const { atoms, box, metadata } = sceneData;

  // Calculate camera distance from box size
  const maxDim = Math.max(box?.a || 10, box?.b || 10, box?.c || 10);

  return (
    <>
      <ambientLight intensity={0.6} />
      <directionalLight position={[maxDim, maxDim, maxDim]} intensity={0.8} />
      <directionalLight position={[-maxDim, -maxDim / 2, maxDim]} intensity={0.3} />
      <AtomInstances atoms={atoms} />
      <OrbitControls
        enableDamping
        dampingFactor={0.1}
        minDistance={2}
        maxDistance={maxDim * 4}
      />
    </>
  );
}

const SceneViewer = ({ sceneData, title }) => {
  if (!sceneData || !sceneData.atoms || sceneData.atoms.length === 0) {
    return (
      <div className="molecule-card">
        <div className="molecule-title">{title || 'Material'}</div>
        <div className="viewer" style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', color: '#666' }}>
          No scene data available
        </div>
      </div>
    );
  }

  const maxDim = Math.max(
    sceneData.box?.a || 10,
    sceneData.box?.b || 10,
    sceneData.box?.c || 10
  );

  return (
    <div className="molecule-card">
      <div className="molecule-title">
        {title && (
          <a
            href={`https://en.wikipedia.org/wiki/${encodeURIComponent(title)}`}
            target="_blank"
            rel="noopener noreferrer"
          >
            {title}
          </a>
        )}
        <span style={{ fontSize: '10px', color: '#888', marginLeft: '8px' }}>
          {sceneData.metadata?.visualization_mode} | {sceneData.metadata?.atom_count} atoms
        </span>
      </div>
      <div className="viewer">
        <Canvas
          camera={{
            position: [maxDim * 0.8, maxDim * 0.6, maxDim * 0.8],
            fov: 50,
            near: 0.1,
            far: maxDim * 10
          }}
          style={{ background: '#000000' }}
        >
          <SceneContent sceneData={sceneData} />
        </Canvas>
      </div>
    </div>
  );
};

export default SceneViewer;
