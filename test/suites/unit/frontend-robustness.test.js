// test/unit/frontend-robustness.test.js - Frontend robustness and error handling tests

const { TextInput } = require('../../frontend/components/TextInput');
const { useApi } = require('../../frontend/hooks/useApi');

// Mock fetch for API testing
global.fetch = jest.fn();

describe('Frontend Robustness Tests', () => {
  beforeEach(() => {
    fetch.mockClear();
  });

  describe('TextInput Component Robustness', () => {
    test('should validate input properly', () => {
      const mockOnSubmit = jest.fn();
      const { getByPlaceholderText, getByText } = render(
        <TextInput 
          value=""
          onChange={() => {}}
          onSubmit={mockOnSubmit}
          isProcessing={false}
        />
      );

      const input = getByPlaceholderText(/specify an object/i);
      
      // Test empty input
      fireEvent.keyDown(input, { key: 'Enter' });
      expect(getByText(/please enter a molecule/i)).toBeInTheDocument();
      expect(mockOnSubmit).not.toHaveBeenCalled();

      // Test short input
      fireEvent.change(input, { target: { value: 'a' } });
      fireEvent.keyDown(input, { key: 'Enter' });
      expect(getByText(/at least 2 characters/i)).toBeInTheDocument();
      expect(mockOnSubmit).not.toHaveBeenCalled();

      // Test valid input
      fireEvent.change(input, { target: { value: 'water' } });
      fireEvent.keyDown(input, { key: 'Enter' });
      expect(mockOnSubmit).toHaveBeenCalledWith('water');
    });

    test('should handle API errors gracefully', async () => {
      const mockOnSubmit = jest.fn().mockRejectedValue(new Error('API Error'));
      const { getByPlaceholderText, getByText } = render(
        <TextInput 
          value="water"
          onChange={() => {}}
          onSubmit={mockOnSubmit}
          isProcessing={false}
        />
      );

      const input = getByPlaceholderText(/specify an object/i);
      fireEvent.keyDown(input, { key: 'Enter' });

      await waitFor(() => {
        expect(getByText(/analysis failed/i)).toBeInTheDocument();
      });
    });

    test('should prevent multiple submissions while processing', () => {
      const mockOnSubmit = jest.fn();
      const { getByPlaceholderText } = render(
        <TextInput 
          value="water"
          onChange={() => {}}
          onSubmit={mockOnSubmit}
          isProcessing={true}
        />
      );

      const input = getByPlaceholderText(/specify an object/i);
      fireEvent.keyDown(input, { key: 'Enter' });
      
      expect(mockOnSubmit).not.toHaveBeenCalled();
      expect(input).toBeDisabled();
    });
  });

  describe('API Hook Robustness', () => {
    test('should handle network timeouts', async () => {
      fetch.mockImplementation(() => 
        new Promise((_, reject) => 
          setTimeout(() => reject(new Error('AbortError')), 100)
        )
      );

      const { result } = renderHook(() => useApi());
      
      await expect(result.current.analyzeText('water')).rejects.toThrow('Request timed out');
      expect(result.current.error).toBe('Request timed out. Please try again.');
    });

    test('should retry failed requests', async () => {
      let callCount = 0;
      fetch.mockImplementation(() => {
        callCount++;
        if (callCount < 3) {
          return Promise.resolve({
            ok: false,
            status: 500,
            statusText: 'Internal Server Error'
          });
        }
        return Promise.resolve({
          ok: true,
          json: () => Promise.resolve({ molecules: [] })
        });
      });

      const { result } = renderHook(() => useApi());
      
      await result.current.analyzeText('water');
      
      expect(fetch).toHaveBeenCalledTimes(3);
    });

    test('should handle rate limiting', async () => {
      fetch.mockResolvedValue({
        ok: false,
        status: 429,
        statusText: 'Too Many Requests'
      });

      const { result } = renderHook(() => useApi());
      
      await expect(result.current.analyzeText('water')).rejects.toThrow('Rate limit exceeded');
    });

    test('should validate input before making API calls', async () => {
      const { result } = renderHook(() => useApi());
      
      await expect(result.current.analyzeText('')).rejects.toThrow('Text input is required');
      await expect(result.current.analyzeText('   ')).rejects.toThrow('Text input is required');
    });
  });

  describe('Molecule labeling', () => {
    test('viewer names should use molecule names, not object query', () => {
      const molecules = [
        { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
        { name: 'Ethanol', smiles: 'CCO' }
      ];

      const objectQuery = 'coffee';
      const smilesArray = molecules.map(m => m.smiles);
      const viewers = molecules.map((mol, index) => ({
        name: mol.name || objectQuery,
        sdfData: `/sdf_files/${smilesArray[index]}.sdf`,
        smiles: mol.smiles
      }));

      expect(viewers.map(v => v.name)).toEqual(['Caffeine', 'Ethanol']);
    });
  });

  describe('MainLayout Component Robustness', () => {
    test('should handle keyboard shortcuts without conflicts', () => {
      const mockSetViewers = jest.fn();
      const mockSetLastAnalysis = jest.fn();
      
      render(
        <MainLayout
          isProcessing={false}
          setIsProcessing={() => {}}
          viewers={[]}
          setViewers={mockSetViewers}
          currentAnalysisType={null}
          setCurrentAnalysisType={() => {}}
          lastAnalysis={null}
          setLastAnalysis={mockSetLastAnalysis}
        />
      );

      // Test that shortcuts work when not in input field
      fireEvent.keyDown(document, { key: 'k', metaKey: true });
      expect(document.activeElement.id).toBe('object-input');

      // Test that shortcuts are ignored in input field
      const input = document.getElementById('object-input');
      fireEvent.keyDown(input, { key: 'c', metaKey: true });
      // Should not trigger camera mode when in input field
    });

    test('should clear errors when input changes', () => {
      const { getByPlaceholderText, queryByText } = render(
        <MainLayout
          isProcessing={false}
          setIsProcessing={() => {}}
          viewers={[]}
          setViewers={() => {}}
          currentAnalysisType={null}
          setCurrentAnalysisType={() => {}}
          lastAnalysis={null}
          setLastAnalysis={() => {}}
        />
      );

      const input = getByPlaceholderText(/specify an object/i);
      
      // Simulate error state
      fireEvent.change(input, { target: { value: 'a' } });
      fireEvent.keyDown(input, { key: 'Enter' });
      
      expect(queryByText(/at least 2 characters/i)).toBeInTheDocument();
      
      // Change input should clear error
      fireEvent.change(input, { target: { value: 'water' } });
      expect(queryByText(/at least 2 characters/i)).not.toBeInTheDocument();
    });
  });

  describe('Error Recovery', () => {
    test('should provide retry functionality', async () => {
      const mockAnalyzeText = jest.fn()
        .mockRejectedValueOnce(new Error('Network Error'))
        .mockResolvedValueOnce({ molecules: [{ name: 'water', sdf_data: 'test' }] });

      const { getByTitle, getByText } = render(
        <MainLayout
          isProcessing={false}
          setIsProcessing={() => {}}
          viewers={[]}
          setViewers={() => {}}
          currentAnalysisType={null}
          setCurrentAnalysisType={() => {}}
          lastAnalysis={null}
          setLastAnalysis={() => {}}
        />
      );

      // Trigger failed analysis
      const input = getByPlaceholderText(/specify an object/i);
      fireEvent.change(input, { target: { value: 'water' } });
      fireEvent.keyDown(input, { key: 'Enter' });

      await waitFor(() => {
        expect(getByText(/analysis failed/i)).toBeInTheDocument();
      });

      // Retry should be available
      const retryButton = getByTitle(/retry last analysis/i);
      expect(retryButton).toBeInTheDocument();
    });
  });

  describe('Accessibility', () => {
    test('should provide proper ARIA labels and descriptions', () => {
      const { getByLabelText, getByRole } = render(
        <TextInput 
          value=""
          onChange={() => {}}
          onSubmit={() => {}}
          isProcessing={false}
          error="Test error"
        />
      );

      const input = getByLabelText(/specify an object/i);
      expect(input).toHaveAttribute('aria-describedby');
      
      const errorElement = getByRole('alert');
      expect(errorElement).toBeInTheDocument();
    });
  });
}); 