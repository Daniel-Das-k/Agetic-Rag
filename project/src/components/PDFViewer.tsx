
import React, { useState, useEffect } from 'react';
import { Document, Page, pdfjs } from 'react-pdf';
import { Loader2 } from 'lucide-react';

// Set up PDF.js worker
pdfjs.GlobalWorkerOptions.workerSrc = `//cdnjs.cloudflare.com/ajax/libs/pdf.js/${pdfjs.version}/pdf.worker.min.js`;

interface PDFViewerProps {
  file: File;
}

export const PDFViewer: React.FC<PDFViewerProps> = ({ file }) => {
  const [numPages, setNumPages] = useState<number | null>(null);
  const [pdfUrl, setPdfUrl] = useState<string | null>(null);
  const [pageWidth, setPageWidth] = useState(800);

  useEffect(() => {
    // Create URL for the uploaded file
    const url = URL.createObjectURL(file);
    setPdfUrl(url);

    // Calculate initial page width
    const containerWidth = Math.min(window.innerWidth - 48, 800);
    setPageWidth(containerWidth);

    // Cleanup URL on component unmount
    return () => {
      URL.revokeObjectURL(url);
    };
  }, [file]);

  useEffect(() => {
    // Handle window resize
    const handleResize = () => {
      const containerWidth = Math.min(window.innerWidth - 48, 800);
      setPageWidth(containerWidth);
    };

    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, []);

  const onDocumentLoadSuccess = ({ numPages }: { numPages: number }) => {
    setNumPages(numPages);
  };

  if (!pdfUrl) {
    return (
      <div className="flex items-center justify-center h-32">
        <Loader2 className="h-8 w-8 animate-spin text-gray-400" />
      </div>
    );
  }

  return (
    <div className="flex flex-col items-center py-6 px-4 max-h-[80vh] overflow-y-auto">
      <Document
        file={pdfUrl}
        onLoadSuccess={onDocumentLoadSuccess}
        loading={
          <div className="flex items-center justify-center h-32">
            <Loader2 className="h-8 w-8 animate-spin text-gray-400" />
          </div>
        }
        error={
          <div className="flex items-center justify-center h-32 text-red-500">
            Error loading PDF. Please try again.
          </div>
        }
      >
        {Array.from(new Array(numPages), (_, index) => (
          <div key={`page_${index + 1}`} className="mb-6">
            <Page
              pageNumber={index + 1}
              renderTextLayer={false}
              className="shadow-xl rounded-lg overflow-hidden"
              width={pageWidth}
              loading={
                <div className="flex items-center justify-center h-32">
                  <Loader2 className="h-8 w-8 animate-spin text-gray-400" />
                </div>
              }
            />
          </div>
        ))}
      </Document>
      {numPages && (
        <div className="sticky bottom-4 bg-white/80 backdrop-blur-sm py-2 px-6 rounded-full shadow-lg border border-gray-100">
          <p className="text-sm font-medium text-gray-600">
            {numPages} page{numPages === 1 ? '' : 's'}
          </p>
        </div>
      )}
    </div>
  );
};