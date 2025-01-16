import chromadb
import os
import autogen
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

def termination_msg(x):
    return isinstance(x, dict) and "TERMINATE" == str(x.get("content", ""))[-9:].upper()

Final_paper = AssistantAgent(
    name="Question_format",
    system_message="""Organize the following questions into a finalized question paper. The paper should include:
1. *Short-Answer Questions (2 Marks)* focusing on Remember and Understand.
2. *Long-Answer Questions (10 Marks)* focusing on Evaluate and Create.
3. *Diagram & Scenario-Based Questions (10 Marks)* focusing on Apply and Create.
Ensure the questions are categorized, formatted appropriately, and distributed according to Bloom’s Taxonomy, with a logical structure.
""",
    llm_config=llm_config
)

short_answer = autogen.AssistantAgent(
    name="2_marks",
    is_termination_msg=termination_msg,
    system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 2-mark short-answer questions focusing on *Remembering* (recalling facts) and *Understanding* (explaining concepts). These questions should test knowledge recall and comprehension, ensuring answers are brief (2-3 sentences).",
    llm_config=llm_config,
)

long_answer = autogen.AssistantAgent(
    name="10_marks",
    is_termination_msg=termination_msg,
    system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 10-mark long-answer questions that require *Evaluation* (critically assessing ideas or arguments) or *Creation* (designing new ideas, developing models). Ensure these questions are multi-part and require a detailed response.",
    llm_config=llm_config,
)

diag_answer = autogen.AssistantAgent(
    name="diag_10_marks",
    is_termination_msg=termination_msg,
    system_message="""Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 10-mark diagram-based and scenario-based questions that:
1. Require students to *apply* their knowledge by drawing or labeling diagrams.
2. Involve real-world scenarios where students *create* solutions or demonstrate how theoretical concepts work in practice.
""",
    llm_config=llm_config,
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
    # system_message="Assistant who has extra content retrieval power for solving difficult problems.",
    code_execution_config=False,
    retrieve_config={
        "model": gemini_config_list[0]["model"],
        "task": "qa",
        "update_context": True,
        "n_results": 3,
        "docs_path": [
            "lesson.pdf"
            # os.path.join(os.path.abspath(""), "..", "website", "docs"),
        ],
        "client": chromadb.PersistentClient(path=CHROMA_DB_PATH),
        "get_or_create": True,
        "overwrite": False,
        "vector_db": vector_db,
        "collection_name": CHROMA_COLLECTION,
        "embedding_function": custom_embedding,  
    },
)

PROBLEM="From the provided document, extract the most relevant topics and subtopics. For each topic, retrieve the corresponding explanations, definitions, or related content. Focus on identifying both *factual* and *conceptual* content that could span across different cognitive levels (e.g., factual recall, application, analysis)."

groupchat = autogen.GroupChat(
    agents=[ragproxyagent, short_answer, long_answer, diag_answer,Final_paper], messages=[], max_round=9, speaker_selection_method="round_robin"
)
manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=llm_config)
ragproxyagent.initiate_chat(
    manager,
    problem=PROBLEM,
    message=ragproxyagent.message_generator,
    n_results=3,
)
