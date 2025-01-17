import PyPDF2

def extract_text_from_pdf(pdf_path):
    try:
        with open(pdf_path, 'rb') as pdf_file:
            pdf_reader = PyPDF2.PdfReader(pdf_file)
            extracted_text = ""
            for page in pdf_reader.pages:
                extracted_text += page.extract_text()
            return extracted_text
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Example Usage
pdf_path = "chap_2.pdf"  # Replace with the path to your PDF file
text = extract_text_from_pdf(pdf_path)

if text:
    print("Extracted Text:")
    print(text)
else:
    print("Failed to extract text.")