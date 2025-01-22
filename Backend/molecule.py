import os
import sys
import base64
from rdkit import Chem
from rdkit.Chem import Draw
from autogen import AssistantAgent, UserProxyAgent
from autogen.agentchat.contrib.multimodal_conversable_agent import MultimodalConversableAgent
from autogen import config_list_from_json
from dotenv import load_dotenv

load_dotenv()

def encode_image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode("utf-8")

# Get all image paths
image_folder = "images"
image_paths = [os.path.join(image_folder, img) for img in os.listdir(image_folder) if img.endswith('.png')]
encoded_images = {os.path.basename(path): encode_image_to_base64(path) for path in image_paths}

# Prepare configuration
gemini_config_list = config_list_from_json(
    "OAI_CONFIG_LIST.json",
    filter_dict={"model": [os.getenv("MODEL")]},
)

llm_config = {
    "config_list": gemini_config_list,
    "cache_seed": None,
    "temperature": 0,
    "timeout": 300,
}

def select_most_educational_image(images):
    """
    Simple function to select the most educational image based on predefined criteria
    This can be replaced with more sophisticated selection logic
    """
    # For this example, we'll choose Paclitaxel as it's a complex molecule with medical significance
    return "Paclitaxel.png"

def main():
    # Directly select the most educational image
    selected_image = select_most_educational_image(list(encoded_images.keys()))
    
    # Print the selected image
    print(f"Selected image: {selected_image}")
    
    # Optional: You can add more processing here if needed
    return selected_image

if __name__ == "__main__":
    main()
