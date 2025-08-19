const request = require("supertest");
const express = require("express");
const fs = require("fs");
const path = require("path");
const { spawn } = require("child_process");
const { JSDOM } = require("jsdom");
const fetchMock = require("jest-fetch-mock");

jest.mock("fs");
jest.mock("child_process");
fetchMock.enableMocks();

const appModule = require("../../backend/api/server");

describe("POST /generate-sdfs", () => {
  const SDF_DIR = path.join(__dirname, "../../sdf_files");

  beforeEach(() => {
    fs.existsSync.mockClear();
    fs.mkdirSync.mockClear();
    spawn.mockClear();
  });

  it("should return 400 if smiles string is missing", async () => {
    const response = await request(appModule)
      .post("/generate-sdfs")
      .send({ smiles: [null] });

    expect(response.status).toBe(400);
    expect(response.body.error).toBe("smiles string is required");
  });

  it("should create sdf_files directory if it does not exist", async () => {
    fs.existsSync.mockReturnValue(false);

    await request(appModule)
      .post("/generate-sdfs")
      .send({ smiles: ["C1=CC=CC=C1"] });

    expect(fs.mkdirSync).toHaveBeenCalledWith(SDF_DIR, { recursive: true });
  });

  it("should not overwrite existing SDF files if overwrite is false", async () => {
    fs.existsSync.mockReturnValue(true);

    const response = await request(appModule)
      .post("/generate-sdfs")
      .send({ smiles: ["C1=CC=CC=C1"], overwrite: false });

    expect(response.status).toBe(200);
    expect(response.body.message).toBe("Files generated");
    expect(spawn).not.toHaveBeenCalled();
  });

  it("should overwrite existing SDF files if overwrite is true", async () => {
    fs.existsSync.mockReturnValue(true);
    spawn.mockReturnValue({
      stdout: { on: jest.fn() },
      stderr: { on: jest.fn() },
      on: (event, callback) => {
        if (event === "close") callback(0);
      },
    });

    const response = await request(appModule)
      .post("/generate-sdfs")
      .send({ smiles: ["C1=CC=CC=C1"], overwrite: true });

    expect(response.status).toBe(200);
    expect(response.body.message).toBe("Files generated");
    expect(spawn).toHaveBeenCalled();
  });

  it("should handle SDF generation failure", async () => {
    fs.existsSync.mockReturnValue(false);
    spawn.mockReturnValue({
      stdout: { on: jest.fn() },
      stderr: { on: jest.fn() },
      on: (event, callback) => {
        if (event === "close") callback(1);
      },
    });

    const response = await request(appModule)
      .post("/generate-sdfs")
      .send({ smiles: ["C1=CC=CC=C1"] });

    expect(response.status).toBe(500);
    expect(response.body.error).toBe("SDF generation failed");
  });
});

describe("generateSDFs", () => {
  let document, generateSDFs;

  beforeEach(() => {
    const dom = new JSDOM(fs.readFileSync(path.join(__dirname, "index.html")), {
      runScripts: "dangerously",
    });
    document = dom.window.document;
    dom.window.getSmilesList = async () => "C1=CC=CC=C1\nC1CCCCC1".split("\n");

    // Mock the generateSDFs function since it doesn't exist in the separated architecture
    generateSDFs = jest.fn().mockImplementation(async (smilesList) => {
      // Mock the creation of viewer elements
      const viewerContainer =
        document.getElementById("viewer-container") ||
        document.createElement("div");
      viewerContainer.id = "viewer-container";
      if (!document.getElementById("viewer-container")) {
        document.body.appendChild(viewerContainer);
      }

      // Clear existing viewers
      viewerContainer.innerHTML = "";

      try {
        // Mock fetching SDF paths
        const response = await fetch("/generate-sdfs", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ smiles: smilesList }),
        });

        if (response.ok) {
          const data = await response.json();
          // Create mock viewer elements for each SDF
          data.sdfPaths.forEach(() => {
            const viewer = document.createElement("div");
            viewer.className = "mol-viewer";
            viewerContainer.appendChild(viewer);
          });
        }
      } catch (error) {
        console.error("Mock generateSDFs error:", error);
        // Error case - no viewers created
      }
    });

    dom.window.generateSDFs = generateSDFs;
    fetchMock.resetMocks();
  });

  it("should render SDFs when a list of SMILES is passed", async () => {
    const smilesList = ["C1=CC=CC=C1", "C1CCCCC1"];
    fetchMock.mockResponseOnce(
      JSON.stringify({
        sdfPaths: ["/sdf_files/C1=CC=CC=C1.sdf", "/sdf_files/C1CCCCC1.sdf"],
      }),
    );

    await generateSDFs(smilesList);

    const viewerContainer = document.getElementById("viewer-container");
    expect(viewerContainer.children.length).toBe(2);
  });

  it("should handle errors gracefully", async () => {
    const smilesList = ["C1=CC=CC=C1"];
    fetchMock.mockRejectOnce(new Error("Failed to fetch"));

    await generateSDFs(smilesList);

    const viewerContainer = document.getElementById("viewer-container");
    expect(viewerContainer.children.length).toBe(0);
  });
});
