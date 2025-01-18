
// import React, { useState, useCallback } from 'react';
// import { pdfjs } from 'react-pdf';
// import { Upload, FileUp, Loader2, GraduationCap } from 'lucide-react';
// import { PDFViewer } from './components/PDFViewer';
// import { Button } from './components/Button';
// import axios from 'axios';

// pdfjs.GlobalWorkerOptions.workerSrc = `//cdnjs.cloudflare.com/ajax/libs/pdf.js/${pdfjs.version}/pdf.worker.min.js`;

// function App() {
//   const [file, setFile] = useState<File | null>(null);
//   const [isDragging, setIsDragging] = useState(false);
//   const [isLoading, setIsLoading] = useState(false);
//   const [analysisData, setAnalysisData] = useState<any | null>(null); 

//   const handleDragOver = useCallback((e: React.DragEvent) => {
//     e.preventDefault();
//     setIsDragging(true);
//   }, []);

//   const handleDragLeave = useCallback((e: React.DragEvent) => {
//     e.preventDefault();
//     setIsDragging(false);
//   }, []);

//   const handleDrop = useCallback((e: React.DragEvent) => {
//     e.preventDefault();
//     setIsDragging(false);
    
//     const droppedFile = e.dataTransfer.files[0];
//     if (droppedFile?.type === 'application/pdf') {
//       setFile(droppedFile);
//     }
//   }, []);

//   const handleFileChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
//     const selectedFile = e.target.files?.[0];
//     if (selectedFile?.type === 'application/pdf') {
//       setFile(selectedFile);
//     }
//   }, []);

//   const handleGenerate = async () => {
//     if (!file) return;
  
//     setIsLoading(true);
//     const formData = new FormData();
//     formData.append('file', file);
  
//     try {
//       const response = await axios.post('http://127.0.0.1:9090/run', formData);
  
//       if (response.status !== 200) {
//         throw new Error('Failed to process PDF');
//       }
  
//       const data = response.data;
//       setAnalysisData(data); 
//       setFile(null); 
//       console.log('Success:', data);
//     } catch (error) {
//       console.error('Error:', error);
//     } finally {
//       setIsLoading(false);
//     }
//   };

//   return (
//     <div className="min-h-screen bg-gradient-to-br from-indigo-50 via-white to-purple-50">
//       {/* Header */}
//       <header className="bg-white/80 backdrop-blur-sm border-b border-indigo-100 sticky top-0 z-50">
//         <div className="max-w-7xl mx-auto px-4 py-4 sm:px-6 lg:px-8">
//           <div className="flex items-center justify-between">
//             <div className="flex items-center space-x-3">
//               <div className="p-2 bg-indigo-600 rounded-lg">
//                 <GraduationCap className="h-6 w-6 text-white" />
//               </div>
//               <div>
//                 <h1 className="text-xl font-semibold text-gray-900">EduSage</h1>
//                 <p className="text-sm text-gray-500">PDF Analysis Tool</p>
//               </div>
//             </div>
//             {file && (
//               <p className="text-sm text-gray-600 bg-gray-100 px-3 py-1 rounded-full">
//                 {file.name}
//               </p>
//             )}
//           </div>
//         </div>
//       </header>

//       <main className="max-w-7xl mx-auto px-4 py-8 sm:px-6 lg:px-8">
//         {/* Show PDF preview or backend content */}
//         {!analysisData ? (
//           !file ? (
//             <div
//               className={`relative h-[70vh] border-2 border-dashed rounded-2xl transition-all duration-200 ${
//                 isDragging
//                   ? 'border-indigo-500 bg-indigo-50/50 scale-[0.99]'
//                   : 'border-gray-300 bg-white/50 hover:border-indigo-300 hover:bg-white/80'
//               }`}
//               onDragOver={handleDragOver}
//               onDragLeave={handleDragLeave}
//               onDrop={handleDrop}
//             >
//               <div className="absolute inset-0 flex flex-col items-center justify-center">
//                 <div className="p-4 bg-white rounded-full shadow-md mb-6">
//                   <Upload className="h-8 w-8 text-indigo-600" />
//                 </div>
//                 <h2 className="text-xl font-medium text-gray-700 mb-2">
//                   Upload your PDF document
//                 </h2>
//                 <p className="text-sm text-gray-500 mb-4">
//                   Drag and drop your file here, or{' '}
//                   <label className="text-indigo-600 hover:text-indigo-500 cursor-pointer font-medium">
//                     browse
//                     <input
//                       type="file"
//                       className="hidden"
//                       accept=".pdf"
//                       onChange={handleFileChange}
//                     />
//                   </label>
//                 </p>
//                 <p className="text-xs text-gray-400">Supported format: PDF</p>
//               </div>
//             </div>
//           ) : (
//             <div className="bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100">
//               <PDFViewer file={file} />
//             </div>
//           )
//         ) : (
//           <div className="bg-white/80 backdrop-blur-sm rounded-2xl shadow-xl border border-gray-100 p-6">
//             <h2 className="text-xl font-medium text-gray-900 mb-4">Analysis Output</h2>
//             <pre className="text-sm text-gray-700">{JSON.stringify(analysisData, null, 2)}</pre>
//           </div>
//         )}
//       </main>

//       {/* Floating Action Button */}
//       {file && (
//         <div className="fixed bottom-8 right-8">
//           <Button
//             onClick={handleGenerate}
//             disabled={isLoading}
//             className="shadow-lg px-6 py-3 text-base"
//           >
//             {isLoading ? (
//               <>
//                 <Loader2 className="h-5 w-5 animate-spin" />
//                 <span className="ml-2">Processing...</span>
//               </>
//             ) : (
//               <>
//                 <FileUp className="h-5 w-5" />
//                 <span className="ml-2">Generate Analysis</span>
//               </>
//             )}
//           </Button>
//         </div>
//       )}
//     </div>
//   );
// }

// export default App;

import React, { useState, useCallback } from 'react';
import { Upload, FileUp, Loader2, GraduationCap, ScrollText } from 'lucide-react';
import axios from 'axios';

interface ExamData {
  title: string;
  instructions: string[];
  total_marks: number;
  sections: {
    section_title: string;
    total_marks: number;
    background?: string;
    supporting_data?: string[];
    questions: {
      question_number: number;
      question_text: string;
      marks: number;
    }[];
  }[];
}

const App = () => {
  const [file, setFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [analysisData, setAnalysisData] = useState<ExamData | null>(null);

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
      setAnalysisData(response.data);
      console.log('Success:', response.data);
      setFile(null);
    } catch (error) {
      console.error('Error:', error);
    } finally {
      setIsLoading(false);
    }
  };

  const renderExamPaper = (data: ExamData) => {
    return (
      <div className="max-w-4xl mx-auto bg-white rounded-xl shadow-lg overflow-hidden">
        {/* Header */}
        <div className="bg-indigo-600 text-white px-6 py-4">
          <div className="flex items-center space-x-3 mb-2">
            <ScrollText className="h-6 w-6" />
            <h1 className="text-2xl font-bold">{data.title}</h1>
          </div>
          <div className="bg-indigo-500/50 rounded-lg p-4 mt-4">
            <h2 className="font-semibold mb-2">Instructions:</h2>
            <ul className="list-disc list-inside space-y-1">
              {data.instructions.map((instruction, index) => (
                <li key={index} className="text-sm">{instruction}</li>
              ))}
            </ul>
          </div>
        </div>

        {/* Sections */}
        <div className="px-6 py-8 space-y-8">
          {data.sections.map((section, sectionIndex) => (
            <div key={sectionIndex} className="space-y-4">
              <div className="flex justify-between items-center border-b border-indigo-100 pb-2">
                <h2 className="text-xl font-semibold text-gray-900">
                  {section.section_title}
                </h2>
                <span className="text-sm text-indigo-600 font-medium">
                  {section.total_marks} marks
                </span>
              </div>

              {/* Case Study Background */}
              {section.background && (
                <div className="bg-gray-50 rounded-lg p-4 mb-4">
                  <p className="text-gray-700">{section.background}</p>
                  {section.supporting_data && (
                    <div className="mt-3">
                      <h4 className="font-medium text-gray-900 mb-2">Supporting Data:</h4>
                      <ul className="list-disc list-inside text-sm text-gray-600">
                        {section.supporting_data.map((data, index) => (
                          <li key={index}>{data}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              )}

              {/* Questions */}
              <div className="space-y-4">
                {section.questions.map((question, questionIndex) => (
                  <div 
                    key={questionIndex}
                    className="p-4 rounded-lg border border-gray-200 hover:border-indigo-200 transition-colors"
                  >
                    <div className="flex justify-between items-start mb-2">
                      <span className="font-medium text-gray-900">
                        Q{question.question_number}.
                      </span>
                      <span className="text-sm text-indigo-600 bg-indigo-50 px-2 py-1 rounded">
                        {question.marks} marks
                      </span>
                    </div>
                    <p className="text-gray-700">{question.question_text}</p>
                  </div>
                ))}
              </div>
            </div>
          ))}
        </div>

        {/* Footer */}
        <div className="bg-gray-50 px-6 py-4 border-t border-gray-200">
          <div className="flex justify-between items-center">
            <span className="text-gray-600">Total Marks:</span>
            <span className="text-lg font-semibold text-indigo-600">{data.total_marks}</span>
          </div>
        </div>
      </div>
    );
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
        {!analysisData ? (
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
            {renderExamPaper(analysisData)}
          </div>
        )}
      </main>

      {/* Floating Action Button */}
      {file && (
        <div className="fixed bottom-8 right-8">
          <button
            onClick={handleGenerate}
            disabled={isLoading}
            className="flex items-center space-x-2 bg-indigo-600 hover:bg-indigo-700 text-white px-6 py-3 rounded-lg shadow-lg disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {isLoading ? (
              <>
                <Loader2 className="h-5 w-5 animate-spin" />
                <span>Processing...</span>
              </>
            ) : (
              <>
                <FileUp className="h-5 w-5" />
                <span>Generate Analysis</span>
              </>
            )}
          </button>
        </div>
      )}
    </div>
  );
};

export default App;