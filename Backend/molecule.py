from rdkit import Chem
from rdkit.Chem import Draw
from autogen import AssistantAgent
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from autogen.agentchat.contrib.m_conversable_agent_multimodal import MultimodalConversableAgent
from autogen import config_list_from_json
import autogen
from autogen import UserProxyAgent
import os
from dotenv import load_dotenv
load_dotenv()

# smiles = "C1=CC=C(C=C1)O"
# mol = Chem.MolFromSmiles(smiles)


# Draw.MolToFile(mol, "Phenol.png")
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


user_proxy = UserProxyAgent(
        name = "user_proxy",
        system_message = """
            You are the UserProxyAgent, responsible for coordinating the retrieval and analysis of images for educational purposes. 
            Your goal is to retrieve images from the folder(images/), analyze them to determine which image best suits the purpose of 
            posing challenging questions to students, and save the selected image to the output folder. 
            Ensure clear communication between the user and the assistant agents, maintaining efficiency throughout the process.
            The Images are in images/ folder location
        """,
         human_input_mode="NEVER",
        llm_config=llm_config
    )

image_anlayser = AssistantAgent(
            name = "image_analyser",
            system_message = """
You are the ImageAnalyserAgent, responsible for analyzing images in the specified 'images/' folder. 
Your task is to evaluate all images in the folder(images/) and select the one best suited for educational purposes. 
The chosen image should be challenging, thought-provoking, and capable of forming the basis for questions that test 
students' skills effectively. Consider factors like visual complexity, uniqueness, and relevance to educational scenarios. 
Once the analysis is complete, provide the path or filename of the selected image to the next designated agent(image_retrieve) for saving.
""",

 human_input_mode="NEVER",
            llm_config=llm_config
        )

retreive_images = AssistantAgent(
            name = "image_retrieve",
            system_message = """
You are the ImageRetrieveAgent, responsible for handling and saving the image selected by the ImageAnalyserAgent. 
You will receive the path or filename of the chosen image from the 'images' folder. Your task is to retrieve this image 
and save it in the local directory at the specified location: 'chosen/1.png'. Ensure that the directory 'chosen' exists, 
creating it if necessary, before saving the image. Confirm the successful completion of the task and provide the saved 
image path as output.
"ALWAYS save the figure in `chosen/result.jpg` file. Tell other agents it is in the <img chosen/result.jpg> file location."
"""

,
            llm_config=llm_config,
            human_input_mode="NEVER"
        )


groupchat = autogen.GroupChat(agents=[user_proxy, image_anlayser, retreive_images], messages=[], max_round=12)
manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=llm_config)
user_proxy.initiate_chat(
    manager, message="Choose the best image from the folder(images/) for posing questions to students."
)