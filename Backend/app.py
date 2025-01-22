import os
import json
from main import QuestionPaperGenerator
from image_q import MolecularCaseStudyGenerator
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

class IntegratedQuestionGenerator:
    def __init__(self, file_path):
        """
        Initialize the integrated question generator with both standard and image-based capabilities.
        
        Args:
            file_path: Path to the input file for question generation
        """
        self.file_path = file_path
        
        # Ensure environment variables are set
        if not os.getenv("MODEL") or not os.getenv("URL") or not os.getenv("EMBEDDING"):
            raise ValueError("""
                Please ensure all required environment variables are set in .env file:
                MODEL = "gemini-2.0-flash-exp"
                URL = "your-google-api-key"
                EMBEDDING = "models/embedding-001"
            """)
            
        # Initialize generators
        print(f"Initializing question generator with file: {file_path}")
        self.question_generator = QuestionPaperGenerator(file_path)
        self.molecular_generator = MolecularCaseStudyGenerator()
        
    def generate_short_answers(self):
        """Generate short answer questions"""
        print("Generating short answer questions...")
        response = self.question_generator.ragproxyagent.initiate_chat(
            self.question_generator.short_answer,
            message="Generate short answer questions."
        )
        return response.chat_history[-1]['content'] if response.chat_history else None

    def generate_long_answers(self):
        """Generate long answer questions"""
        print("Generating long answer questions...")
        response = self.question_generator.ragproxyagent.initiate_chat(
            self.question_generator.long_answer,
            message="Generate long answer questions."
        )
        return response.chat_history[-1]['content'] if response.chat_history else None

    def generate_case_study(self):
        """Generate case study based on availability of images"""
        images_folder = os.path.join(os.path.dirname(self.file_path), "images")
        print(f"Checking for images folder at: {images_folder}")
        
        try:
            if os.path.exists(images_folder) and os.path.isdir(images_folder):
                # Use image-based case study
                image_files = [f for f in os.listdir(images_folder) if f.endswith(('.png', '.jpg', '.jpeg'))]
                if image_files:
                    print(f"Found images: {image_files}")
                    # Using specific image for case study
                    image_path = "images/Benzene.png"
                    print(f"Using image: {image_path}")
                    case_study = self.molecular_generator.generate_case_study(image_path)
                    # Ensure case study is properly formatted
                    if isinstance(case_study, str):
                        try:
                            # Try to parse as JSON if it's a string
                            case_study = json.loads(case_study)
                        except json.JSONDecodeError:
                            pass
                    return case_study
                else:
                    print("No images found in images folder, falling back to standard case study")
                    return self.generate_standard_case_study()
            else:
                print("No images folder found, using standard case study")
                return self.generate_standard_case_study()
        except Exception as e:
            print(f"Error in case study generation: {str(e)}")
            return self.generate_standard_case_study()

    def generate_standard_case_study(self):
        """Generate a standard (non-image-based) case study"""
        print("Generating standard case study...")
        response = self.question_generator.ragproxyagent.initiate_chat(
            self.question_generator.case_study,
            message="Generate a case study question."
        )
        return response.chat_history[-1]['content'] if response.chat_history else None

    def combine_and_format(self, short_answers, long_answers, case_study):
        """Combine all sections and format as JSON"""
        print("Combining all sections...")
        
        # Ensure case_study is a string
        if isinstance(case_study, dict):
            case_study = json.dumps(case_study)
            
        combined = self.question_generator.ragproxyagent.initiate_chat(
            self.question_generator.final_paper,
            message=f"""Combine these sections:
            Short Answer: {short_answers}
            Long Answer: {long_answers}
            Case Study: {case_study}"""
        )
        
        if not combined.chat_history:
            raise Exception("Failed to combine sections")
            
        print("Formatting final paper as JSON...")
        json_formatted = self.question_generator.ragproxyagent.initiate_chat(
            self.question_generator.json_formater,
            message=combined.chat_history[-1]['content']
        )
        
        if not json_formatted.chat_history:
            raise Exception("Failed to format JSON")
            
        final_content = json_formatted.chat_history[-1]['content']
        
        # Ensure the output is valid JSON
        try:
            if isinstance(final_content, str):
                return json.loads(final_content)
            return final_content
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON: {str(e)}")
            return {"error": "Failed to generate valid JSON output"}

    def generate_questions(self):
        """
        Generate a complete question paper with short answers, long answers, and case studies.
        Returns a properly formatted JSON response.
        """
        try:
            # Generate each section independently
            short_answers = self.generate_short_answers()
            if not short_answers:
                raise Exception("Failed to generate short answers")
                
            long_answers = self.generate_long_answers()
            if not long_answers:
                raise Exception("Failed to generate long answers")
                
            case_study = self.generate_case_study()
            if not case_study:
                raise Exception("Failed to generate case study")

            # Combine and format all sections
            final_paper = self.combine_and_format(short_answers, long_answers, case_study)
            return final_paper
            
        except Exception as e:
            print(f"Error generating questions: {str(e)}")
            return {"error": str(e)}

def main():
    # Specify your input file path
    file_path = "Electric_charges_and_fields.pdf"
    
    try:
        generator = IntegratedQuestionGenerator(file_path)
        question_paper = generator.generate_questions()
        
        # Ensure we have a valid response
        if isinstance(question_paper, dict) and "error" in question_paper:
            print(f"\nError in question paper generation: {question_paper['error']}")
            return None
            
        print("\nGenerated Question Paper:")
        print(json.dumps(question_paper, indent=2))
        return question_paper
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return None

if __name__ == "__main__":
    main()
