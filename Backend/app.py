import os
from dotenv import load_dotenv
from flask import Flask, request, jsonify
from flask_cors import CORS
from qa import QuestionPaperGenerator
from image_case import ImageCaseStudyGenerator

app = Flask(__name__)
CORS(app)

load_dotenv()

required_env_vars = ["MODEL", "URL", "EMBEDDING"]
missing_vars = [var for var in required_env_vars if not os.getenv(var)]
if missing_vars:
    raise EnvironmentError(f"Missing required environment variables: {', '.join(missing_vars)}")

def process_file(file_path):
    image_dir = "images/"
    print("Image directory:", image_dir)
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

    if os.path.exists(image_dir) and os.path.isdir(image_dir):
        image_files = [f for f in os.listdir(image_dir) if f.lower().endswith(('.png', '.jpg', '.jpeg'))]
        if image_files:
            image_path = os.path.join(image_dir, image_files[0])
            print(f"Using image file: {image_path}")
            image_generator = ImageCaseStudyGenerator(file_path, image_path)
            return image_generator.generate_paper(image_path)
    
    print("No images found, generating text-based questions")
    text_generator = QuestionPaperGenerator(file_path)
    return text_generator.generate_paper()


@app.route('/generate', methods=['POST'])
def generate_questions():
    try:
        if 'file' not in request.files:
            return jsonify({"error": "No file provided"}), 400
            
        file = request.files['file']
        if not file.filename:
            return jsonify({"error": "No file selected"}), 400
            
        print("Processing file:", file.filename)

        upload_dir = "/tmp/uploads"
        os.makedirs(upload_dir, exist_ok=True)

        file_path = os.path.join(upload_dir, file.filename)
        file.save(file_path)
        print("File saved to:", file_path)
        
        result = process_file(file_path)

        import shutil
        shutil.rmtree(upload_dir)
        
        return jsonify({
            "final_output": result
        })
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return jsonify({
            "error": str(e),
            "message": "Please check your input files and try again"
        }), 500

if __name__ == "__main__":
    app.run(debug=True, port=9090)