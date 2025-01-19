
import React, { useState, useCallback } from 'react';
import { pdfjs } from 'react-pdf';
import { Upload, FileUp, Loader2, GraduationCap } from 'lucide-react';
import { PDFViewer } from './components/PDFViewer';
import { Button } from './components/ui/Button';
pdfjs.GlobalWorkerOptions.workerSrc = `//cdnjs.cloudflare.com/ajax/libs/pdf.js/${pdfjs.version}/pdf.worker.min.js`;

function App() {
  const [file, setFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [analysisData, setAnalysisData] = useState<any | null>(null);
  const [showGeneratedPDF, setShowGeneratedPDF] = useState(false);
  const [error, setError] = useState<string | null>(null);

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
      setError(null);
    }
  }, []);

  const handleFileChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const selectedFile = e.target.files?.[0];
    if (selectedFile?.type === 'application/pdf') {
      setFile(selectedFile);
      setError(null);
    }
  }, []);

  const handleGenerate = async () => {
    if (!file) return;
  
    setIsLoading(true);
    setError(null);
  
    try {
      
      const formData = new FormData();
      formData.append('file', file);
  
    
      const response = await fetch('http://localhost:9090/run', {
        method: 'POST',
        body: formData,
      });
  
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
  
      const data = await response.json();
  
     
      let datas = data.final_output;
  
      if (datas && typeof datas === 'string') {
        datas = datas.trim().slice(7, -3); 
      }
  
      let parsedData;
      try {
        parsedData = JSON.parse(datas);
      } catch (parseError) {
        console.error('JSON Parsing Error:', parseError);
        throw new Error('Failed to parse the server response.');
      }
  
      console.log(typeof parsedData); 
  
      setAnalysisData(parsedData);
      setShowGeneratedPDF(true);
    } catch (error) {
      console.error('Error:', error);
      setError('Failed to analyze PDF. Please try again.');
  

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
        {error && (
          <div className="mb-4 p-4 bg-red-50 border border-red-200 rounded-lg">
            <p className="text-red-600">{error}</p>
          </div>
        )}
        
        <div className="relative flex overflow-hidden">
          {/* PDF Viewer Section */}
          <div
            className={`flex-1 transition-all duration-500 ease-in-out transform 
              ${showGeneratedPDF ? 'w-1/3 -translate-x-2/3 opacity-50 scale-95' : 'w-full translate-x-0 opacity-100'}
              bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100 p-6`}
          >
            {!file ? (
              <div
                className={`relative h-[70vh] border-2 border-dashed rounded-2xl transition-all duration-200 
                  ${isDragging ? 'border-indigo-500 bg-indigo-50/50 scale-[0.99]' : 'border-gray-300 bg-white/50 hover:border-indigo-300 hover:bg-white/80'}`}
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
              <PDFViewer file={file} />
            )}
          </div>

          {/* Generated PDF Section */}
          {analysisData && (
            console.log(analysisData),
            <div
              className={`absolute top-0 right-0 flex-1 transition-all duration-500 ease-in-out transform
                ${showGeneratedPDF ? 'translate-x-0 opacity-100 w-2/3' : 'translate-x-full opacity-0 w-0'}
                bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100 p-6 overflow-auto max-h-[90vh]`}
            >
              <div className="prose max-w-none">
                <h1 className="text-2xl font-bold text-center mb-6">{analysisData.title}</h1>
                
                {/* Instructions */}
                <div className="mb-8">
                  <h2 className="text-lg font-semibold mb-2">Instructions:</h2>
                  <ul className="list-disc pl-5">
                    {analysisData.instructions.map((instruction: string, index: number) => (
                      <li key={index} className="text-gray-700">{instruction}</li>
                    ))}
                  </ul>
                </div>
                
                {/* Sections */}
                {analysisData.sections.map((section: any, sectionIndex: number) => (
                  <div key={sectionIndex} className="mb-8">
                    <div className="flex justify-between items-center mb-4">
                      <h2 className="text-xl font-semibold">{section.section_title}</h2>
                      <span className="text-gray-600">Total Marks: {section.total_marks}</span>
                    </div>
                    
                    {/* Case Study Background */}
                    {section.background && (
                      <div className="bg-gray-50 p-4 rounded-lg mb-4">
                        <p className="text-gray-700">{section.background}</p>
                        {section.supporting_data && (
                          <ul className="list-disc pl-5 mt-2">
                            {section.supporting_data.map((data: string, index: number) => (
                              <li key={index} className="text-gray-600">{data}</li>
                            ))}
                          </ul>
                        )}
                      </div>
                    )}
                    
                    {/* Questions */}
                    <div className="space-y-4">
                      {section.questions.map((question: any, questionIndex: number) => (
                        <div key={questionIndex} className="flex">
                          <div className="mr-4 font-medium min-w-[2rem]">
                            {question.question_number}.
                          </div>
                          <div className="flex-grow">
                            <div className="flex justify-between">
                              <p className="text-gray-800">{question.question_text}</p>
                              <span className="text-gray-500 ml-4">[{question.marks} marks]</span>
                            </div>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
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
