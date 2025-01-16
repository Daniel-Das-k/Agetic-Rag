import chromadb
import os
from autogen import AssistantAgent
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from autogen import config_list_from_json
from langchain_google_genai import GoogleGenerativeAIEmbeddings
from langchain_chroma import Chroma

class CustomEmbeddingFunction:
    def __init__(self, embedding_model):
        self.embedding_model = embedding_model

    def __call__(self, input: list) -> list:
        """
        Convert input texts into embeddings using the embedding model.

        Args:
            input (list): List of strings to embed.

        Returns:
            list: List of embeddings, where each embedding is a list of floats.
        """
        return self.embedding_model.embed_documents(input)

gemini_config_list = config_list_from_json(
    "OAI_CONFIG_LIST.json",
    filter_dict={"model": ["gemini-2.0-flash-exp"]},
)

llm_config = {
    "config_list": gemini_config_list,
    "seed": 53,
    "temperature": 0,
    "timeout": 300,
}

assistant = AssistantAgent(
    name="assistant",
    system_message="You are a helpful assistant.",
    llm_config=llm_config
)

CHROMA_DB_PATH = "/tmp/chromadb"
CHROMA_COLLECTION = "autogen-docs-test"

chroma_client = chromadb.PersistentClient(path=CHROMA_DB_PATH)
collection = chroma_client.get_or_create_collection(name=CHROMA_COLLECTION)

embeddings = GoogleGenerativeAIEmbeddings(
    model="models/embedding-001",
    google_api_key="AIzaSyBt3JTEUjrNFb392jP-2YExCCSOlv8ev7A"
)
custom_embedding = CustomEmbeddingFunction(embedding_model=embeddings)

vector_db = Chroma(embedding_function=custom_embedding)

ragproxyagent = RetrieveUserProxyAgent(
    name="ragproxyagent",
    human_input_mode="NEVER",
    llm_config=llm_config,
    code_execution_config=False,
    retrieve_config={
        "model": gemini_config_list[0]["model"],
        "task": "qa",
        "update_context": True,
        "n_results": 3,
        "docs_path": [
            "Kenesis.pdf",
            os.path.join(os.path.abspath(""), "..", "website", "docs"),
        ],
        "client": chromadb.PersistentClient(path=CHROMA_DB_PATH),
        "get_or_create": True,
        "overwrite": False,
        "vector_db": vector_db,
        "collection_name": CHROMA_COLLECTION,
        "embedding_function": custom_embedding,  
    },
)

code_problem = (
    "Generate 10 questions according to BLOOM taxonomy from the provided document")

ragproxyagent.initiate_chat(
    assistant, message=ragproxyagent.message_generator, problem=code_problem)
