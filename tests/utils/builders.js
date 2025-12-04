/**
 * Test data builders using the builder pattern
 * Provides flexible and consistent test data creation
 */

/**
 * Base builder class
 */
class Builder {
  constructor() {
    this.data = {};
  }
  
  with(key, value) {
    this.data[key] = value;
    return this;
  }
  
  build() {
    return { ...this.data };
  }
}

/**
 * Molecule data builder
 */
class MoleculeBuilder extends Builder {
  constructor() {
    super();
    // Set defaults
    this.data = {
      id: Math.random().toString(36).substr(2, 9),
      name: 'Ethanol',
      smiles: 'CCO',
      formula: 'C2H6O',
      mass: 46.07,
      createdAt: new Date().toISOString()
    };
  }
  
  withName(name) {
    this.data.name = name;
    return this;
  }
  
  withSmiles(smiles) {
    this.data.smiles = smiles;
    return this;
  }
  
  withFormula(formula) {
    this.data.formula = formula;
    return this;
  }
  
  withMass(mass) {
    this.data.mass = mass;
    return this;
  }
  
  withProperties(properties) {
    this.data.properties = properties;
    return this;
  }
  
  asWater() {
    this.data.name = 'Water';
    this.data.smiles = 'O';
    this.data.formula = 'H2O';
    this.data.mass = 18.015;
    return this;
  }
  
  asCaffeine() {
    this.data.name = 'Caffeine';
    this.data.smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C';
    this.data.formula = 'C8H10N4O2';
    this.data.mass = 194.19;
    return this;
  }
  
  asInvalid() {
    this.data.name = '';
    this.data.smiles = 'INVALID_SMILES_@#$%';
    return this;
  }
}

/**
 * API request builder
 */
class RequestBuilder extends Builder {
  constructor() {
    super();
    this.data = {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json'
      }
    };
  }
  
  withMethod(method) {
    this.data.method = method;
    return this;
  }
  
  withBody(body) {
    this.data.body = typeof body === 'object' ? JSON.stringify(body) : body;
    return this;
  }
  
  withHeaders(headers) {
    this.data.headers = { ...this.data.headers, ...headers };
    return this;
  }
  
  withAuth(token) {
    this.data.headers.Authorization = `Bearer ${token}`;
    return this;
  }
  
  asPost() {
    this.data.method = 'POST';
    return this;
  }
  
  asPut() {
    this.data.method = 'PUT';
    return this;
  }
  
  asDelete() {
    this.data.method = 'DELETE';
    return this;
  }
  
  asMultipart() {
    delete this.data.headers['Content-Type'];
    return this;
  }
}

/**
 * API response builder
 */
class ResponseBuilder extends Builder {
  constructor() {
    super();
    this.data = {
      status: 200,
      success: true,
      data: null,
      message: 'Success'
    };
  }
  
  withStatus(status) {
    this.data.status = status;
    this.data.success = status >= 200 && status < 300;
    return this;
  }
  
  withData(data) {
    this.data.data = data;
    return this;
  }
  
  withMessage(message) {
    this.data.message = message;
    return this;
  }
  
  withError(error) {
    this.data.success = false;
    this.data.error = error;
    return this;
  }
  
  asError(status = 500, message = 'Internal Server Error') {
    this.data.status = status;
    this.data.success = false;
    this.data.message = message;
    delete this.data.data;
    return this;
  }
  
  asNotFound() {
    return this.asError(404, 'Not Found');
  }
  
  asUnauthorized() {
    return this.asError(401, 'Unauthorized');
  }
  
  asBadRequest(errors = {}) {
    this.asError(400, 'Bad Request');
    this.data.errors = errors;
    return this;
  }
}

/**
 * User/Session builder
 */
class UserBuilder extends Builder {
  constructor() {
    super();
    this.data = {
      id: Math.random().toString(36).substr(2, 9),
      email: `test-${Date.now()}@example.com`,
      name: 'Test User',
      role: 'user',
      createdAt: new Date().toISOString(),
      isActive: true
    };
  }
  
  withEmail(email) {
    this.data.email = email;
    return this;
  }
  
  withName(name) {
    this.data.name = name;
    return this;
  }
  
  asAdmin() {
    this.data.role = 'admin';
    return this;
  }
  
  asGuest() {
    this.data.role = 'guest';
    this.data.email = null;
    return this;
  }
  
  asInactive() {
    this.data.isActive = false;
    return this;
  }
  
  withApiKey(key) {
    this.data.apiKey = key;
    return this;
  }
}

/**
 * File upload builder
 */
class FileBuilder extends Builder {
  constructor() {
    super();
    this.data = {
      filename: 'test.jpg',
      mimetype: 'image/jpeg',
      size: 1024,
      buffer: Buffer.from('fake-image-data')
    };
  }
  
  withName(filename) {
    this.data.filename = filename;
    return this;
  }
  
  withMimetype(mimetype) {
    this.data.mimetype = mimetype;
    return this;
  }
  
  withSize(size) {
    this.data.size = size;
    return this;
  }
  
  withBuffer(buffer) {
    this.data.buffer = buffer;
    return this;
  }
  
  asPng() {
    this.data.filename = 'test.png';
    this.data.mimetype = 'image/png';
    return this;
  }
  
  asPdf() {
    this.data.filename = 'test.pdf';
    this.data.mimetype = 'application/pdf';
    return this;
  }
  
  asSdf() {
    this.data.filename = 'molecule.sdf';
    this.data.mimetype = 'chemical/x-mdl-sdfile';
    this.data.buffer = Buffer.from(`
$$$$
Ethanol
  -OEChem-01012100003D

  9  8  0     0  0  0  0  0  0999 V2000
    1.2660    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3660    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
    `);
    return this;
  }
  
  asLargeFile() {
    this.data.size = 10 * 1024 * 1024; // 10MB
    this.data.buffer = Buffer.alloc(this.data.size);
    return this;
  }
}

/**
 * Test scenario builder for complex test setups
 */
class ScenarioBuilder {
  constructor() {
    this.molecules = [];
    this.users = [];
    this.requests = [];
    this.responses = [];
  }
  
  withMolecules(count, builderFn) {
    for (let i = 0; i < count; i++) {
      const builder = new MoleculeBuilder();
      if (builderFn) builderFn(builder, i);
      this.molecules.push(builder.build());
    }
    return this;
  }
  
  withUsers(count, builderFn) {
    for (let i = 0; i < count; i++) {
      const builder = new UserBuilder();
      if (builderFn) builderFn(builder, i);
      this.users.push(builder.build());
    }
    return this;
  }
  
  withRequests(count, builderFn) {
    for (let i = 0; i < count; i++) {
      const builder = new RequestBuilder();
      if (builderFn) builderFn(builder, i);
      this.requests.push(builder.build());
    }
    return this;
  }
  
  build() {
    return {
      molecules: this.molecules,
      users: this.users,
      requests: this.requests,
      responses: this.responses
    };
  }
}

/**
 * Factory functions for quick creation
 */
const create = {
  molecule: () => new MoleculeBuilder(),
  request: () => new RequestBuilder(),
  response: () => new ResponseBuilder(),
  user: () => new UserBuilder(),
  file: () => new FileBuilder(),
  scenario: () => new ScenarioBuilder()
};

module.exports = {
  MoleculeBuilder,
  RequestBuilder,
  ResponseBuilder,
  UserBuilder,
  FileBuilder,
  ScenarioBuilder,
  create
};
