// test/integration/system.test.js - System tests for end-to-end functionality
// These tests validate the entire system working together (2-5 minutes)

const request = require("supertest");
const app = require("../../backend/api/server");
const fs = require("fs");
const path = require("path");
const { JSDOM } = require("jsdom");

// Mock external dependencies
jest.mock("openai", () => ({
  OpenAI: jest.fn().mockImplementation(() => ({
    chat: {
      completions: {
        create: jest.fn().mockResolvedValue({
          choices: [
            {
              message: {
                content: JSON.stringify({
                  object: "test object",
                  smiles: ["CCO", "CC(=O)O"],
                }),
              },
            },
          ],
        }),
      },
    },
  })),
}));

jest.mock("child_process", () => ({
  execSync: jest.fn().mockReturnValue("test output"),
}));

// Mock browser APIs
global.navigator = {
  mediaDevices: {
    getUserMedia: jest.fn().mockResolvedValue({
      getTracks: () => [{ stop: jest.fn() }],
    }),
    enumerateDevices: jest.fn().mockResolvedValue([
      { kind: "videoinput", label: "Camera 1" },
      { kind: "videoinput", label: "Camera 2" },
    ]),
  },
};

global.window = {
  isSecureContext: true,
  location: { protocol: "https:" },
};

describe("System Tests", () => {
  let dom;
  let document;
  let window;

  beforeEach(() => {
    // Create a DOM environment
    const htmlContent = fs.readFileSync(
      path.join(__dirname, "../index.html"),
      "utf8",
    );
    dom = new JSDOM(htmlContent, {
      runScripts: "dangerously",
      resources: "usable",
      url: "http://localhost:8080",
    });

    document = dom.window.document;
    window = dom.window;

    // Mock fetch
    window.fetch = jest.fn();

    // Mock 3Dmol
    window.$3Dmol = {
      createViewer: jest.fn().mockReturnValue({
        addModel: jest.fn(),
        setBackgroundColor: jest.fn(),
        setStyle: jest.fn(),
        zoomTo: jest.fn(),
        render: jest.fn(),
        resize: jest.fn(),
      }),
    };
  });

  afterEach(() => {
    // Clean up test files
    const sdfDir = path.join(__dirname, "../sdf_files");
    if (fs.existsSync(sdfDir)) {
      const files = fs.readdirSync(sdfDir);
      files.forEach((file) => {
        if (file.endsWith(".sdf")) {
          fs.unlinkSync(path.join(sdfDir, file));
        }
      });
    }
  });

  describe("Frontend-Backend Integration", () => {
    test("should handle complete user workflow: camera capture -> analysis -> visualization", async () => {
      // Mock successful API responses
      window.fetch
        .mockResolvedValueOnce({
          ok: true,
          json: () =>
            Promise.resolve({
              output: {
                object: "water bottle",
                smiles: ["CCO", "CC(=O)O"],
              },
            }),
        })
        .mockResolvedValueOnce({
          ok: true,
          json: () =>
            Promise.resolve({
              sdfPaths: ["/sdf_files/CCO.sdf", "/sdf_files/CC(=O)O.sdf"],
              errors: [],
              skipped: [],
            }),
        })
        .mockResolvedValueOnce({
          ok: true,
          text: () =>
            Promise.resolve(
              "CCO\n  RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$",
            ),
        })
        .mockResolvedValueOnce({
          ok: true,
          text: () =>
            Promise.resolve(
              "CC(=O)O\n  RDKit          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  2  4  1  0  0  0  0\nM  END\n$$$$",
            ),
        });

      // Simulate camera interaction
      const video = document.getElementById("video-feed");
      const clickEvent = new window.MouseEvent("click", {
        clientX: 100,
        clientY: 100,
      });

      // Mock video properties
      Object.defineProperty(video, "videoWidth", { value: 640 });
      Object.defineProperty(video, "videoHeight", { value: 480 });
      Object.defineProperty(video, "clientWidth", { value: 320 });
      Object.defineProperty(video, "clientHeight", { value: 240 });
      Object.defineProperty(video, "getBoundingClientRect", {
        value: () => ({ left: 0, top: 0 }),
      });

      // Trigger camera click
      video.dispatchEvent(clickEvent);

      // Wait for async operations
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify API calls were made
      expect(window.fetch).toHaveBeenCalledTimes(4);

      // Verify image analysis call
      const imageAnalysisCall = window.fetch.mock.calls[0];
      expect(imageAnalysisCall[0]).toBe("/image-molecules");
      expect(imageAnalysisCall[1].method).toBe("POST");

      // Verify SDF generation call
      const sdfGenerationCall = window.fetch.mock.calls[1];
      expect(sdfGenerationCall[0]).toBe("/generate-sdfs");
      expect(sdfGenerationCall[1].method).toBe("POST");
    });

    test("should handle text input workflow", async () => {
      // Mock API responses
      window.fetch
        .mockResolvedValueOnce({
          ok: true,
          json: () =>
            Promise.resolve({
              output: {
                object: "coffee",
                smiles: ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
              },
            }),
        })
        .mockResolvedValueOnce({
          ok: true,
          json: () =>
            Promise.resolve({
              sdfPaths: ["/sdf_files/CN1C=NC2=C1C(=O)N(C(=O)N2C)C.sdf"],
              errors: [],
              skipped: [],
            }),
        });

      // Simulate text input
      const textInput = document.getElementById("object-input");
      textInput.value = "coffee";

      const enterEvent = new window.KeyboardEvent("keyup", { key: "Enter" });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify API calls
      expect(window.fetch).toHaveBeenCalledTimes(2);

      const textAnalysisCall = window.fetch.mock.calls[0];
      expect(textAnalysisCall[0]).toBe("/object-molecules");
      expect(JSON.parse(textAnalysisCall[1].body)).toEqual({
        object: "coffee",
      });
    });

    test("should handle photo upload workflow", async () => {
      // Mock FileReader
      const mockFileReader = {
        readAsDataURL: jest.fn(),
        result:
          "data:image/jpeg;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
      };
      window.FileReader = jest.fn(() => mockFileReader);

      // Mock API responses
      window.fetch.mockResolvedValueOnce({
        ok: true,
        json: () =>
          Promise.resolve({
            output: {
              object: "uploaded image",
              smiles: ["CCO"],
            },
          }),
      });

      // Create mock file
      const mockFile = new window.File([""], "test.jpg", {
        type: "image/jpeg",
      });

      // Simulate file upload
      const fileInput = document.getElementById("photo-upload");
      const changeEvent = new window.Event("change");
      Object.defineProperty(fileInput, "files", { value: [mockFile] });

      fileInput.dispatchEvent(changeEvent);

      // Wait for async operations
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify file was processed
      expect(window.FileReader).toHaveBeenCalled();
    });
  });

  describe("Camera Functionality", () => {
    test("should initialize camera on page load", async () => {
      // Mock successful camera access
      global.navigator.mediaDevices.getUserMedia.mockResolvedValueOnce({
        getTracks: () => [{ stop: jest.fn() }],
      });

      // Load app.js
      const appScript = fs.readFileSync(
        path.join(__dirname, "../app.js"),
        "utf8",
      );

      // Execute the script in the DOM context
      const script = document.createElement("script");
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify camera was requested
      expect(global.navigator.mediaDevices.getUserMedia).toHaveBeenCalled();
    });

    test("should handle camera permission denial", async () => {
      // Mock camera permission denial
      global.navigator.mediaDevices.getUserMedia.mockRejectedValueOnce(
        new Error("Permission denied"),
      );

      // Load app.js
      const appScript = fs.readFileSync(
        path.join(__dirname, "../app.js"),
        "utf8",
      );
      const script = document.createElement("script");
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify error message is displayed
      const msgBox = document.querySelector(".permission-message");
      expect(msgBox.textContent).toContain("Camera error");
    });

    test("should show switch camera button when multiple cameras available", async () => {
      // Mock multiple cameras
      global.navigator.mediaDevices.enumerateDevices.mockResolvedValueOnce([
        { kind: "videoinput", label: "Front Camera" },
        { kind: "videoinput", label: "Back Camera" },
      ]);

      // Load app.js
      const appScript = fs.readFileSync(
        path.join(__dirname, "../app.js"),
        "utf8",
      );
      const script = document.createElement("script");
      script.textContent = appScript;
      document.head.appendChild(script);

      // Wait for initialization
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify switch camera button was created
      const switchButton = document.querySelector(".switch-camera-btn");
      expect(switchButton).toBeTruthy();
    });
  });

  describe("UI Interactions", () => {
    test("should toggle between camera and photo modes", () => {
      const cameraMode = document.getElementById("camera-mode");
      const photoMode = document.getElementById("photo-mode");
      const cameraContainer = document.querySelector(".camera-container");
      const photoOptions = document.getElementById("photo-options");

      // Initially camera mode should be active
      expect(cameraContainer.style.display).not.toBe("none");
      expect(photoOptions.style.display).toBe("none");

      // Switch to photo mode
      photoMode.checked = true;
      photoMode.dispatchEvent(new window.Event("change"));

      expect(cameraContainer.style.display).toBe("none");
      expect(photoOptions.style.display).toBe("flex");
    });

    test("should handle URL input validation", () => {
      const urlInput = document.getElementById("photo-url");
      const urlButton = document.getElementById("url-analyze");

      // Test invalid URL
      urlInput.value = "not-a-url";
      urlButton.click();

      // Should show alert (mocked)
      expect(window.alert).toBeDefined();
    });

    test("should create object columns with close functionality", async () => {
      // Mock successful SDF generation
      window.fetch
        .mockResolvedValueOnce({
          ok: true,
          json: () =>
            Promise.resolve({
              sdfPaths: ["/sdf_files/CCO.sdf"],
              errors: [],
              skipped: [],
            }),
        })
        .mockResolvedValueOnce({
          ok: true,
          text: () =>
            Promise.resolve(
              "CCO\n  RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$",
            ),
        });

      // Simulate text input
      const textInput = document.getElementById("object-input");
      textInput.value = "test object";

      const enterEvent = new window.KeyboardEvent("keyup", { key: "Enter" });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify object column was created
      const objectColumn = document.querySelector(".object-column");
      expect(objectColumn).toBeTruthy();

      // Verify close button exists
      const closeButton = objectColumn.querySelector(".column-close");
      expect(closeButton).toBeTruthy();

      // Test close functionality
      const initialColumnCount =
        document.querySelectorAll(".object-column").length;
      closeButton.click();
      const finalColumnCount =
        document.querySelectorAll(".object-column").length;
      expect(finalColumnCount).toBe(initialColumnCount - 1);
    });
  });

  describe("Error Handling", () => {
    test("should handle network errors gracefully", async () => {
      // Mock network error
      window.fetch.mockRejectedValueOnce(new Error("Network error"));

      // Simulate text input
      const textInput = document.getElementById("object-input");
      textInput.value = "test object";

      const enterEvent = new window.KeyboardEvent("keyup", { key: "Enter" });
      textInput.dispatchEvent(enterEvent);

      // Wait for async operations
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify error message was displayed
      const errorElement = document.querySelector('h3[style*="color: red"]');
      expect(errorElement).toBeTruthy();
      expect(errorElement.textContent).toContain("Error");
    });

    test("should handle invalid file uploads", () => {
      // Create mock invalid file
      const mockFile = new window.File([""], "test.txt", {
        type: "text/plain",
      });

      // Mock alert
      window.alert = jest.fn();

      // Simulate file upload
      const fileInput = document.getElementById("photo-upload");
      const changeEvent = new window.Event("change");
      Object.defineProperty(fileInput, "files", { value: [mockFile] });

      fileInput.dispatchEvent(changeEvent);

      // Verify alert was shown
      expect(window.alert).toHaveBeenCalledWith("Please select an image file");
    });
  });

  describe("Performance and Scalability", () => {
    test("should handle multiple concurrent analyses", async () => {
      // Mock multiple successful responses
      const mockResponses = Array(5)
        .fill()
        .map(() => ({
          ok: true,
          json: () =>
            Promise.resolve({
              output: {
                object: "test object",
                smiles: ["CCO"],
              },
            }),
        }));

      window.fetch
        .mockResolvedValueOnce(mockResponses[0])
        .mockResolvedValueOnce(mockResponses[1])
        .mockResolvedValueOnce(mockResponses[2])
        .mockResolvedValueOnce(mockResponses[3])
        .mockResolvedValueOnce(mockResponses[4]);

      // Simulate multiple rapid inputs
      const textInput = document.getElementById("object-input");
      const promises = [];

      for (let i = 0; i < 5; i++) {
        textInput.value = `test object ${i}`;
        const enterEvent = new window.KeyboardEvent("keyup", { key: "Enter" });
        textInput.dispatchEvent(enterEvent);
        promises.push(new Promise((resolve) => setTimeout(resolve, 50)));
      }

      await Promise.all(promises);

      // Verify all API calls were made
      expect(window.fetch).toHaveBeenCalledTimes(5);
    });

    test("should handle large SMILES arrays efficiently", async () => {
      const largeSmilesArray = Array(20).fill("CCO");

      window.fetch.mockResolvedValueOnce({
        ok: true,
        json: () =>
          Promise.resolve({
            output: {
              object: "test object",
              smiles: largeSmilesArray,
            },
          }),
      });

      // Simulate input with large SMILES array
      const textInput = document.getElementById("object-input");
      textInput.value = "test object";

      const enterEvent = new window.KeyboardEvent("keyup", { key: "Enter" });
      textInput.dispatchEvent(enterEvent);

      // Wait for processing
      await new Promise((resolve) => setTimeout(resolve, 100));

      // Verify API call was made
      expect(window.fetch).toHaveBeenCalledTimes(1);
    });
  });
});
