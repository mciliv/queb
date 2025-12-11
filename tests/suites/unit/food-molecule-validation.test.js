// test/unit/food-molecule-validation.test.js - Food-specific molecule validation

const MolecularProcessor = require("../../../src/server/services/molecular-processor");

// Essential molecules found in common foods with their accurate SMILES
const FOOD_MOLECULES = {
  "tomato": {
    essential: [
      { name: "lycopene", smiles: "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C", description: "red pigment antioxidant" },
      { name: "vitamin C", smiles: "C([C@H]([C@H]([C@@H](C(=O)O)O)O)O)O", description: "ascorbic acid" },
      { name: "citric acid", smiles: "C(C(=O)O)C(CC(=O)O)(C(=O)O)O", description: "natural preservative" },
      { name: "beta-carotene", smiles: "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C", description: "vitamin A precursor" }
    ],
    common: [
      { name: "glucose", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", description: "natural sugar" },
      { name: "fructose", smiles: "C([C@@H]([C@@H]([C@H]([C@H](CO)O)O)O)O)O", description: "fruit sugar" }
    ]
  },
  
  "apple": {
    essential: [
      { name: "quercetin", smiles: "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", description: "flavonoid antioxidant" },
      { name: "chlorogenic acid", smiles: "C1=CC(=C(C=C1/C=C/C(=O)O)O)O", description: "polyphenol antioxidant" },
      { name: "malic acid", smiles: "C(C(C(=O)O)O)C(=O)O", description: "gives tart flavor" },
      { name: "pectin (D-galacturonic acid)", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)C(=O)O)O)O)O)O", description: "fiber component" }
    ],
    common: [
      { name: "fructose", smiles: "C([C@@H]([C@@H]([C@H]([C@H](CO)O)O)O)O)O", description: "primary sugar" },
      { name: "glucose", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", description: "secondary sugar" }
    ]
  },
  
  "coffee": {
    essential: [
      { name: "caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", description: "stimulant alkaloid" },
      { name: "chlorogenic acid", smiles: "C1=CC(=C(C=C1/C=C/C(=O)O)O)O", description: "antioxidant" },
      { name: "quinic acid", smiles: "C([C@H]1[C@@H]([C@H]([C@@H]([C@H](C1)O)O)O)O)O", description: "organic acid" },
      { name: "trigonelline", smiles: "C[N+]1=CC=CC=C1C(=O)O", description: "alkaloid precursor" }
    ],
    common: [
      { name: "water", smiles: "O", description: "primary component" }
    ]
  },
  
  "chocolate": {
    essential: [
      { name: "theobromine", smiles: "CN1C=NC2=C1C(=O)NC(=O)N2C", description: "mild stimulant" },
      { name: "caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", description: "stimulant" },
      { name: "epicatechin", smiles: "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O", description: "flavonoid antioxidant" },
      { name: "anandamide", smiles: "CCCCCCCCCCCCCCCCCCCC(=O)NCCO", description: "bliss molecule" }
    ],
    common: [
      { name: "sucrose", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@@H](O2)CO)O)O)CO)O)O)O)O", description: "added sugar" }
    ]
  },
  
  "orange": {
    essential: [
      { name: "vitamin C", smiles: "C([C@H]([C@H]([C@@H](C(=O)O)O)O)O)O", description: "ascorbic acid" },
      { name: "limonene", smiles: "CC1=CCC(CC1)C(=C)C", description: "citrus oil terpene" },
      { name: "hesperidin", smiles: "COC1=CC(=C2C(=C1)OC(=CC2=O)C3=CC(=C(C=C3)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)C)O)O)O)O)O)O)O)O", description: "citrus flavonoid" },
      { name: "citric acid", smiles: "C(C(=O)O)C(CC(=O)O)(C(=O)O)O", description: "natural acid" }
    ],
    common: [
      { name: "fructose", smiles: "C([C@@H]([C@@H]([C@H]([C@H](CO)O)O)O)O)O", description: "natural sugar" },
      { name: "glucose", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", description: "natural sugar" }
    ]
  },
  
  "spinach": {
    essential: [
      { name: "lutein", smiles: "CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC1C(O1)(C)C)C)C=CC=C(C)C=CC=C(C)C=CC2C(=CC(CC2(C)C)O)C", description: "carotenoid for eye health" },
      { name: "folate", smiles: "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NC2=NC3=C(C=NC=N3)N=C2N", description: "vitamin B9" },
      { name: "iron complex", smiles: "[Fe+2]", description: "essential mineral" },
      { name: "vitamin K1", smiles: "CCC[C@H](C)CCC[C@H](C)CCCC(C)CCCC(C)CCCC(C)C", description: "phylloquinone" }
    ],
    common: [
      { name: "water", smiles: "O", description: "primary component" },
      { name: "cellulose unit", smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", description: "fiber" }
    ]
  }
};

describe("Food Molecule Validation Tests", () => {
  let processor;

  beforeEach(() => {
    processor = new MolecularProcessor();
  });

  describe("Essential Food Molecule Recognition", () => {
    test("should validate all essential tomato molecules", () => {
      const tomato = FOOD_MOLECULES.tomato;
      
      tomato.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      tomato.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });

    test("should validate all essential apple molecules", () => {
      const apple = FOOD_MOLECULES.apple;
      
      apple.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      apple.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });

    test("should validate all essential coffee molecules", () => {
      const coffee = FOOD_MOLECULES.coffee;
      
      coffee.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      coffee.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });

    test("should validate all essential chocolate molecules", () => {
      const chocolate = FOOD_MOLECULES.chocolate;
      
      chocolate.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      chocolate.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });

    test("should validate all essential orange molecules", () => {
      const orange = FOOD_MOLECULES.orange;
      
      orange.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      orange.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });

    test("should validate all essential spinach molecules", () => {
      const spinach = FOOD_MOLECULES.spinach;
      
      spinach.essential.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
      
      spinach.common.forEach(molecule => {
        expect(processor.isValidSmiles(molecule.smiles)).toBe(true);
      });
    });
  });

  describe("Food-Specific Molecule Processing", () => {
    test("should process tomato essential molecules without length restrictions", async () => {
      const tomato = FOOD_MOLECULES.tomato;
      const allMolecules = [...tomato.essential, ...tomato.common].map(m => m.smiles);
      
      const result = await processor.processSmiles(allMolecules);
      
      // Should not skip any for length (lycopene and beta-carotene are very long)
      expect(result.skipped.length).toBe(0);
      
      // Lycopene specifically should not be skipped
      const lycopeneSkipped = result.skipped.some(s => s.includes("CC(=CCCC(=CC=CC(=CC=CC=C(C)"));
      expect(lycopeneSkipped).toBe(false);
    });

    test("should process complex food antioxidants", async () => {
      const antioxidants = [
        FOOD_MOLECULES.tomato.essential.find(m => m.name === "lycopene").smiles,
        FOOD_MOLECULES.apple.essential.find(m => m.name === "quercetin").smiles,
        FOOD_MOLECULES.chocolate.essential.find(m => m.name === "epicatechin").smiles,
        FOOD_MOLECULES.orange.essential.find(m => m.name === "hesperidin").smiles,
      ];
      
      const result = await processor.processSmiles(antioxidants);
      
      // Should not skip any for format issues
      expect(result.skipped.length).toBe(0);
      
      // Should attempt to process all
      expect(result.sdfPaths.length + result.errors.length).toBe(antioxidants.length);
    });

    test("should handle vitamins and essential nutrients", async () => {
      const vitamins = [
        FOOD_MOLECULES.tomato.essential.find(m => m.name === "vitamin C").smiles,
        FOOD_MOLECULES.orange.essential.find(m => m.name === "vitamin C").smiles,
        FOOD_MOLECULES.spinach.essential.find(m => m.name === "folate").smiles,
      ];
      
      const result = await processor.processSmiles(vitamins);
      
      // All should be valid format
      expect(result.skipped.length).toBe(0);
    });
  });

  describe("Food Chemistry Accuracy", () => {
    test("should distinguish between similar molecules in different foods", () => {
      // Caffeine in coffee vs theobromine in chocolate
      const caffeine = FOOD_MOLECULES.coffee.essential.find(m => m.name === "caffeine").smiles;
      const theobromine = FOOD_MOLECULES.chocolate.essential.find(m => m.name === "theobromine").smiles;
      
      expect(caffeine).not.toBe(theobromine);
      expect(processor.isValidSmiles(caffeine)).toBe(true);
      expect(processor.isValidSmiles(theobromine)).toBe(true);
    });

    test("should handle stereochemistry in natural sugars", () => {
      const glucose = FOOD_MOLECULES.apple.common.find(m => m.name === "glucose").smiles;
      const fructose = FOOD_MOLECULES.apple.common.find(m => m.name === "fructose").smiles;
      
      // Both should be valid but different
      expect(processor.isValidSmiles(glucose)).toBe(true);
      expect(processor.isValidSmiles(fructose)).toBe(true);
      expect(glucose).not.toBe(fructose);
    });

    test("should validate complex natural product structures", () => {
      const complexMolecules = [
        FOOD_MOLECULES.orange.essential.find(m => m.name === "hesperidin").smiles, // Complex glycoside
        FOOD_MOLECULES.spinach.essential.find(m => m.name === "lutein").smiles, // Complex carotenoid
        FOOD_MOLECULES.chocolate.essential.find(m => m.name === "anandamide").smiles, // Long chain amide
      ];
      
      complexMolecules.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });
  });

  describe("Nutritional Completeness Validation", () => {
    test("should identify all major nutrient classes in foods", () => {
      const nutrients = {
        carbohydrates: [
          FOOD_MOLECULES.apple.common.find(m => m.name === "glucose").smiles,
          FOOD_MOLECULES.apple.common.find(m => m.name === "fructose").smiles,
        ],
        vitamins: [
          FOOD_MOLECULES.tomato.essential.find(m => m.name === "vitamin C").smiles,
          FOOD_MOLECULES.spinach.essential.find(m => m.name === "folate").smiles,
        ],
        antioxidants: [
          FOOD_MOLECULES.tomato.essential.find(m => m.name === "lycopene").smiles,
          FOOD_MOLECULES.apple.essential.find(m => m.name === "quercetin").smiles,
        ],
        alkaloids: [
          FOOD_MOLECULES.coffee.essential.find(m => m.name === "caffeine").smiles,
          FOOD_MOLECULES.chocolate.essential.find(m => m.name === "theobromine").smiles,
        ]
      };
      
      Object.entries(nutrients).forEach(([category, molecules]) => {
        molecules.forEach(smiles => {
          expect(processor.isValidSmiles(smiles)).toBe(true);
        });
      });
    });

    test("should process complete food profiles efficiently", async () => {
      // Test a complete food profile (all molecules from one food)
      const tomatoProfile = [...FOOD_MOLECULES.tomato.essential, ...FOOD_MOLECULES.tomato.common]
        .map(m => m.smiles);
      
      const startTime = Date.now();
      const result = await processor.processSmiles(tomatoProfile);
      const endTime = Date.now();
      
      expect(endTime - startTime).toBeLessThan(15000); // Should complete within 15 seconds
      expect(result.skipped.length).toBe(0); // All should be valid format
    });
  });

  describe("Food Safety and Validation", () => {
    test("should reject obviously non-food molecules", () => {
      const nonFoodMolecules = [
        "H2O", // molecular formula (not SMILES)
        "C6H12O6", // molecular formula
        "NaCl", // salt formula
        "", // empty
        "N/A", // placeholder
      ];
      
      nonFoodMolecules.forEach(molecule => {
        expect(processor.isValidSmiles(molecule)).toBe(false);
      });
    });

    test("should validate natural vs synthetic molecule distinction", () => {
      // Both natural and synthetic vitamin C should have same SMILES
      const naturalVitC = FOOD_MOLECULES.tomato.essential.find(m => m.name === "vitamin C").smiles;
      const syntheticVitC = "C([C@H]([C@H]([C@@H](C(=O)O)O)O)O)O"; // Same structure
      
      expect(naturalVitC).toBe(syntheticVitC);
      expect(processor.isValidSmiles(naturalVitC)).toBe(true);
    });
  });
});

// Export food molecule data for other tests
module.exports = {
  FOOD_MOLECULES,
  getFoodMolecules: (foodName) => FOOD_MOLECULES[foodName],
  getAllFoodSmiles: () => {
    const allSmiles = [];
    Object.values(FOOD_MOLECULES).forEach(food => {
      [...food.essential, ...food.common].forEach(molecule => {
        allSmiles.push(molecule.smiles);
      });
    });
    return allSmiles;
  },
  validateFoodMolecules: async (foodName) => {
    const processor = new MolecularProcessor();
    const food = FOOD_MOLECULES[foodName];
    if (!food) return null;
    
    const allMolecules = [...food.essential, ...food.common].map(m => m.smiles);
    return await processor.processSmiles(allMolecules);
  }
}; 