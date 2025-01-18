
import React, { useState, useCallback } from 'react';
import { pdfjs } from 'react-pdf';
import { Upload, FileUp, Loader2, GraduationCap } from 'lucide-react';
import { PDFViewer } from './components/PDFViewer';
import { Button } from './components/Button';
import axios from 'axios';

pdfjs.GlobalWorkerOptions.workerSrc = `//cdnjs.cloudflare.com/ajax/libs/pdf.js/${pdfjs.version}/pdf.worker.min.js`;

function App() {
  const [file, setFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [analysisData, setAnalysisData] = useState<any | null>(null); 

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
  }, []);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
    
    const droppedFile = e.dataTransfer.files[0];
    if (droppedFile?.type === 'application/pdf') {
      setFile(droppedFile);
    }
  }, []);

  const handleFileChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = e.target.files?.[0];
    if (selectedFile?.type === 'application/pdf') {
      setFile(selectedFile);
    }
  }, []);

  const handleGenerate = async () => {
    if (!file) return;
  
    setIsLoading(true);
    const formData = new FormData();
    formData.append('file', file);
  
    try {
      const response = await axios.post('http://127.0.0.1:9090/run', formData);
  
      if (response.status !== 200) {
        throw new Error('Failed to process PDF');
      }
  
      const data = response.data;
      setAnalysisData(data); 
      setFile(null); 
      console.log('Success:', data);
    } catch (error) {
      console.error('Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-indigo-50 via-white to-purple-50">
      {/* Header */}
      <header className="bg-white/80 backdrop-blur-sm border-b border-indigo-100 sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-4 py-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-3">
              <div className="p-2 bg-indigo-600 rounded-lg">
                <GraduationCap className="h-6 w-6 text-white" />
              </div>
              <div>
                <h1 className="text-xl font-semibold text-gray-900">EduSage</h1>
                <p className="text-sm text-gray-500">PDF Analysis Tool</p>
              </div>
            </div>
            {file && (
              <p className="text-sm text-gray-600 bg-gray-100 px-3 py-1 rounded-full">
                {file.name}
              </p>
            )}
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 py-8 sm:px-6 lg:px-8">
        {/* Show PDF preview or backend content */}
        {!analysisData ? (
          !file ? (
            <div
              className={`relative h-[70vh] border-2 border-dashed rounded-2xl transition-all duration-200 ${
                isDragging
                  ? 'border-indigo-500 bg-indigo-50/50 scale-[0.99]'
                  : 'border-gray-300 bg-white/50 hover:border-indigo-300 hover:bg-white/80'
              }`}
              onDragOver={handleDragOver}
              onDragLeave={handleDragLeave}
              onDrop={handleDrop}
            >
              <div className="absolute inset-0 flex flex-col items-center justify-center">
                <div className="p-4 bg-white rounded-full shadow-md mb-6">
                  <Upload className="h-8 w-8 text-indigo-600" />
                </div>
                <h2 className="text-xl font-medium text-gray-700 mb-2">
                  Upload your PDF document
                </h2>
                <p className="text-sm text-gray-500 mb-4">
                  Drag and drop your file here, or{' '}
                  <label className="text-indigo-600 hover:text-indigo-500 cursor-pointer font-medium">
                    browse
                    <input
                      type="file"
                      className="hidden"
                      accept=".pdf"
                      onChange={handleFileChange}
                    />
                  </label>
                </p>
                <p className="text-xs text-gray-400">Supported format: PDF</p>
              </div>
            </div>
          ) : (
            <div className="bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100">
              <PDFViewer file={file} />
            </div>
          )
        ) : (
          <div className="bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100 p-6">
            <h2 className="text-xl font-medium text-gray-900 mb-4">Analysis Output</h2>
            <pre className="text-sm text-gray-700">{JSON.stringify(analysisData, null, 2)}</pre>
          </div>
        )}
      </main>

      {/* Floating Action Button */}
      {file && (
        <div className="fixed bottom-8 right-8">
          <Button
            onClick={handleGenerate}
            disabled={isLoading}
            className="shadow-lg px-6 py-3 text-base"
          >
            {isLoading ? (
              <>
                <Loader2 className="h-5 w-5 animate-spin" />
                <span className="ml-2">Processing...</span>
              </>
            ) : (
              <>
                <FileUp className="h-5 w-5" />
                <span className="ml-2">Generate Analysis</span>
              </>
            )}
          </Button>
        </div>
      )}
    </div>
  );
}

export default App;
