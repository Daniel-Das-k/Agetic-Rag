import React, { useState } from 'react';
import { Document, Page } from 'react-pdf';

interface PDFViewerProps {
  file: File;
}

export const PDFViewer: React.FC<PDFViewerProps> = ({ file }) => {
  const [numPages, setNumPages] = useState<number | null>(null);

  const onDocumentLoadSuccess = ({ numPages }: { numPages: number }) => {
    setNumPages(numPages);
  };

  return (
    <div className="flex flex-col items-center py-6 px-4 max-h-[80vh] overflow-y-auto scrollbar-thin">
      <Document
        file={file}
        onLoadSuccess={onDocumentLoadSuccess}
        className="max-w-full"
        loading={
          <div className="flex items-center justify-center h-32">
            <div className="animate-pulse text-gray-400">Loading PDF...</div>
          </div>
        }
      >
        {Array.from(new Array(numPages), (_, index) => (
          <div key={`page_${index + 1}`} className="mb-6">
            <Page
              pageNumber={index + 1}
              renderTextLayer={false}
              className="shadow-xl rounded-lg overflow-hidden"
              width={Math.min(window.innerWidth - 48, 800)}
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