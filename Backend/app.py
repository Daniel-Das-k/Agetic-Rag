import os
from dotenv import load_dotenv
from main import QuestionPaperGenerator
from image_q import MolecularCaseStudyGenerator

load_dotenv()

required_env_vars = ["MODEL", "URL", "EMBEDDING"]
missing_vars = [var for var in required_env_vars if not os.getenv(var)]
if missing_vars:
    raise EnvironmentError(f"Missing required environment variables: {', '.join(missing_vars)}")

def process_file(file_path):

    images_dir = os.path.join(os.path.dirname(file_path), "images")
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    if os.path.exists(images_dir) and any(os.path.isfile(os.path.join(images_dir, f)) for f in os.listdir(images_dir)):
        image_files = [f for f in os.listdir(images_dir) if f.lower().endswith(('.png', '.jpg', '.jpeg'))]
        if image_files:
            print(image_files[1],"**"* 1000)
            image_path = os.path.join(images_dir, image_files[1])
            print(f"Using image file: {image_path}")
            image_generator = MolecularCaseStudyGenerator(file_path, image_path)
            return image_generator.generate_paper(image_path)
    
    print("No images found, generating text-based questions")
    text_generator = QuestionPaperGenerator(file_path)
    return text_generator.generate_paper()

if __name__ == "__main__":

    try:
        file_path = "Electric_charges_and_fields.pdf"
        result = process_file(file_path)
        print("\nGenerated Output:")
        print(result)
    except Exception as e:
        print(f"Error: {str(e)}")