import os
import base64
import json
import autogen
from autogen import AssistantAgent, UserProxyAgent
from autogen import config_list_from_json
from dotenv import load_dotenv

class MolecularCaseStudyGenerator:

    def __init__(self, config_path="OAI_CONFIG_LIST.json"):
        # Load environment variables
        load_dotenv()

        # Configuration setup
        self.config_path = config_path
        self.model = os.getenv("MODEL")
        self.llm_config = self._create_llm_config()

    def _create_llm_config(self):
        return {
            "config_list": config_list_from_json(
                self.config_path,
                filter_dict={"model": [self.model]}
            ),
            "cache_seed": 40,
            "temperature": 0,
            "timeout": 180,
        }

    def encode_image_to_base64(self, image_path):

        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode("utf-8")

    def _create_molecular_analyst_agent(self, encoded_image):
        """
        Create the Molecular Structure Analyst agent.
        
        :param encoded_image: Base64 encoded image
        :return: AssistantAgent for molecular structure analysis
        """
        return AssistantAgent(
            name="Molecular_Structure_Analyst",
            system_message=f"""
            You are an advanced molecular structure analysis agent specializing in Bloom's Taxonomy 
            APPLY and CREATE levels of learning.

            Bloom's Taxonomy Guidelines:
            - APPLY Level Objectives:
              * Translate structural observations into practical scenarios
              * Demonstrate how molecular features influence real-world applications
              * Develop problem-solving questions that require applying molecular knowledge

            - CREATE Level Objectives:
              * Design innovative scenarios that extend beyond direct observation
              * Encourage hypothetical thinking and novel applications
              * Generate questions that challenge learners to synthesize and innovate

            Image Analysis Requirements:
            - Carefully analyze the provided molecular structure image: {encoded_image[:100]}...
            
            Specific Analysis Approach:
            1. APPLY: Transform structural insights into practical challenges
               - How can the molecular structure inform industrial or research applications?
               - What real-world problems can be solved using this molecular understanding?

            2. CREATE: Develop forward-thinking, innovative scenarios
               - Propose novel research directions inspired by the molecular structure
               - Design hypothetical applications that push current technological boundaries

            Constraints:
            - Questions must directly derive from the molecular structure image
            - Demonstrate scientific creativity and critical thinking
            - Balance between rigorous analysis and imaginative exploration
            - You can include as many questions as needed(not necessarily)
            """,
            llm_config=self.llm_config
        )

    def _create_json_formatter_agent(self, image_path):
        """
        Create the JSON Formatter agent.
        
        :param image_path: Path to the molecular structure image
        :return: AssistantAgent for JSON formatting
        """
        return AssistantAgent(
            name="JSON_Formatter",
            system_message=f"""
            You are a precise JSON formatting agent with STRICT requirements:

            MANDATORY OUTPUT FORMAT:
            {{
                "section_title": "Molecular Structure: [Specific Molecule] Analysis",
                "total_marks": 15,
                "image": "[molecule_image.png]",
                "questions": [
                    {{
                        "question_number": 1,
                        "question_text": "APPLY-level question with specific scenario, chemical context, and clear application",
                        "marks": X,
                    
                    }},
                    {{
                        "question_number": 2,
                        "question_text": "CREATE-level question requiring innovative thinking and hypothesis generation",
                        "marks": Y,
                      
                    }},
                    ....
                ]
            }}

            STRICT FORMATTING RULES:
            1. TOTAL MARKS MUST EQUAL 15
            2. Include BOTH APPLY and CREATE level questions
            3. Each question MUST have:
               - Clear question text
               - Specific marks allocation
            4. Use the molecule's specific name in section_title
            5. Include the actual image filename

            CRITICAL CONSTRAINTS:
            - NO additional text or explanation
            - PURE JSON output
            - VALIDATE JSON structure before outputting
            - Use the specific molecule image: {os.path.basename(image_path)}
            """,
            llm_config=self.llm_config
        )

    def generate_case_study(self, image_path):

        encoded_image = self.encode_image_to_base64(image_path)

        molecular_analyst = self._create_molecular_analyst_agent(encoded_image)
        json_formatter = self._create_json_formatter_agent(image_path)
        user_proxy = UserProxyAgent(
            name="Scientific_Coordinator",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=2,
            system_message=f"Coordinate precise case study generation for the molecular structure image: {os.path.basename(image_path)}. ENSURE STRICT JSON OUTPUT FORMAT.",
            llm_config=self.llm_config
        )

        # Create Group Chat
        groupchat = autogen.GroupChat(
            agents=[user_proxy, molecular_analyst, json_formatter],
            messages=[],
            max_round=6,
            speaker_selection_method="round_robin"
        )
        manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=self.llm_config)

        # Problem statement
        problem = f"""
        Generate a comprehensive, JSON-structured case study for {os.path.basename(image_path)}
        
        STRICT REQUIREMENTS:
        - Use Bloom's Taxonomy APPLY and CREATE levels
        - Total marks: 15
        - Include both APPLY and CREATE questions
        - Specific marks allocation
        
        Image: {os.path.basename(image_path)}
        Base64 Preview: {encoded_image[:50]}...
        """

        # Initiate the case study generation
        response = user_proxy.initiate_chat(
            manager,
            message=problem,
            n_results=1
        )

        # Extract and validate JSON content
        for message in response.chat_history:
            try:
                # Remove any code block markers
                content = message['content'].replace('```json', '').replace('```', '').strip()
                
                # Parse and validate JSON
                case_study = json.loads(content)
                
                # Additional validation
                assert 'section_title' in case_study
                assert 'total_marks' in case_study and case_study['total_marks'] == 15
                assert 'questions' in case_study
                
                return case_study
            except (json.JSONDecodeError, KeyError, AssertionError) as e:
                print(f"JSON parsing error: {e}")
                continue

        # Fallback if no valid JSON found
        return {"error": "Unable to generate case study"}

def main():
 
    image_path = "images/aceticAcid.png"
    generator = MolecularCaseStudyGenerator()
    case_study = generator.generate_case_study(image_path)
    print(case_study)
    return case_study

if __name__ == "__main__":
    main()
