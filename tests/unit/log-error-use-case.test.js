const LogErrorUseCase = require('../../src/server/services/log-error-use-case');

describe('LogErrorUseCase - Clean Architecture', () => {
  let mockLogErrorPort;
  let logErrorUseCase;

  beforeEach(() => {
    // MOCK PORT: Test the domain logic by mocking infrastructure
    mockLogErrorPort = {
      log: jest.fn()
    };

    logErrorUseCase = new LogErrorUseCase({ logErrorPort: mockLogErrorPort });
  });

  test('should map error type to warn level', async () => {
    // DOMAIN LOGIC TEST: Verify type mapping rules
    await logErrorUseCase.execute({
      type: 'warn',
      message: 'Test warning',
      source: 'frontend'
    });

    expect(mockLogErrorPort.log).toHaveBeenCalledWith({
      level: 'warn',
      message: '[frontend] warn: Test warning',
      metadata: expect.objectContaining({
        type: 'warn',
        message: 'Test warning',
        source: 'frontend'
      })
    });
  });

  test('should map unknown type to info level', async () => {
    // DOMAIN LOGIC TEST: Verify fallback behavior
    await logErrorUseCase.execute({
      type: 'unknown',
      message: 'Unknown type',
      source: 'test'
    });

    expect(mockLogErrorPort.log).toHaveBeenCalledWith({
      level: 'info',
      message: '[test] unknown: Unknown type',
      metadata: expect.any(Object)
    });
  });

  test('should format message with default source', async () => {
    // DOMAIN LOGIC TEST: Verify message formatting rules
    await logErrorUseCase.execute({
      type: 'error',
      message: 'No source provided'
      // source defaults to 'frontend'
    });

    expect(mockLogErrorPort.log).toHaveBeenCalledWith({
      level: 'error',
      message: '[frontend] error: No source provided',
      metadata: expect.any(Object)
    });
  });

  test('should truncate stack traces in metadata', async () => {
    // DOMAIN LOGIC TEST: Verify stack truncation rules
    const longStack = 'Error: something\n'.repeat(100); // 2000+ chars

    await logErrorUseCase.execute({
      type: 'error',
      message: 'Long stack',
      stack: longStack
    });

    const call = mockLogErrorPort.log.mock.calls[0][0];
    expect(call.metadata.stack.length).toBeLessThanOrEqual(500);
    expect(call.metadata.stack.endsWith('...')).toBe(false); // Should be cleanly truncated
  });

  test('should propagate port errors', async () => {
    // INFRASTRUCTURE FAILURE TEST: Verify error propagation
    mockLogErrorPort.log.mockRejectedValue(new Error('Logger failed'));

    await expect(logErrorUseCase.execute({
      type: 'error',
      message: 'Test'
    })).rejects.toThrow('Logger failed');
  });
});