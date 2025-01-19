
# To install required packages:
# pip install pyautogen==0.2.9 panel==1.3.8 google-generativeai
import autogen
import panel as pn
import os
import time
import uuid
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
import asyncio
from langchain_google_genai import GoogleGenerativeAIEmbeddings
from langchain_chroma import Chroma
from autogen import config_list_from_json
from langchain_chroma import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
import chromadb
from autogen import AssistantAgent
import google.generativeai as genai
from dotenv import load_dotenv

load_dotenv()
# Set up Gemini configuration
config_list = config_list_from_json(
    "OAI_CONFIG_LIST.json",
    filter_dict={"model": ["gemini-2.0-flash-exp"]},
)

# Configure the LLM settings
llm_config = {
    "config_list": config_list,
    "cache_seed": None,
    "temperature": 0,
    "timeout": 300,
}

input_future = None

class CustomEmbeddingFunction:
    def __init__(self, embedding_model):
        self.embedding_model = embedding_model

    def __call__(self, input: list) -> list:
        return self.embedding_model.embed_documents(input)
    
unique_id=uuid.uuid4() 
CHROMA_DB_PATH = f"/tmp/chroma.db-{unique_id}"
CHROMA_COLLECTION = "autogen-docs-test"
text_splitter = RecursiveCharacterTextSplitter(separators=["\n\n", "\n", "\r", "\t"])
    
chroma_client = chromadb.PersistentClient(path=CHROMA_DB_PATH)
collection = chroma_client.get_or_create_collection(name=CHROMA_COLLECTION)

embeddings = GoogleGenerativeAIEmbeddings(
    model=os.getenv("EMBEDDING"),
    google_api_key=os.getenv("URL"),
)
custom_embedding = CustomEmbeddingFunction(embedding_model=embeddings)
vector_db = Chroma(embedding_function=custom_embedding)

class MyConversableAgent(autogen.ConversableAgent):
    async def a_get_human_input(self, prompt: str) -> str:
        global input_future
        chat_interface.send(prompt, user="System", respond=False)
        if input_future is None or input_future.done():
            input_future = asyncio.Future()
        await input_future
        input_value = input_future.result()
        input_future = None
        return input_value


def termination_msg(x):
        return isinstance(x, dict) and "TERMINATE" == str(x.get("content", ""))[-9:].upper()

# Update user_proxy setup to include RAG question approval
user_proxy = MyConversableAgent(
    name="Admin",
    is_termination_msg=lambda x: x.get("content", "").rstrip().endswith("exit"),
    system_message="""You are a human admin. Your role is to approve or modify the questions generated by the RAG proxy. Collaborate with other agents to finalize the output.""",
    code_execution_config=False,
    human_input_mode="ALWAYS",
    llm_config=llm_config,
)

# Update RetrieveUserProxyAgent to call user_proxy for approval
ragproxyagent = RetrieveUserProxyAgent(
    name="ragproxyagent",
    human_input_mode="NEVER",  # Enable human intervention
    llm_config=llm_config,
    code_execution_config=False,
    retrieve_config={
        "model": config_list[0]["model"],
        "task": "qa",
        "update_context": True,
        "n_results": 3,
        "docs_path": ["chap_2.pdf"],
        "chunk_token_size": 3000,
        "chunk_mode": "one_line",
        "client": chroma_client,
        "get_or_create": True,
        "overwrite": True,
        "vector_db": vector_db,
        "collection_name": CHROMA_COLLECTION,
        "embedding_function": custom_embedding,
        "custom_text_split_function": text_splitter.split_text,
    },
)




final_paper = AssistantAgent(
            name="Question_format",
            human_input_mode="ALWAYS",
            system_message="""Organize the following questions into a finalized question paper. The paper should include:
            Instuctions: Provide Basic instuctions to the students to write the exam and the marks distribution.
            1. *Short-Answer Questions (2 Marks)* focusing on Remember and Understand.(10 Questions)
            2. *Long-Answer Questions (13 Marks)* focusing on Evaluate and Create.(5 Questions)
            3. *Case Study (15 Marks)* focusing on Apply and Create.(1 Question)
            Ensure the questions are categorized, formatted appropriately, and distributed according to Bloom's Taxonomy, with a logical structure.""",
            llm_config=llm_config
        )

short_answer = autogen.AssistantAgent(
            name="2_marks",
            is_termination_msg=termination_msg,
            human_input_mode="ALWAYS",
            system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 2-mark short-answer questions focusing on *Remembering* and *Understanding*.",
            llm_config=llm_config,
        )
        
long_answer = autogen.AssistantAgent(
            name="10_marks",
            is_termination_msg=termination_msg,
            human_input_mode="ALWAYS",
            system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 13-mark long-answer questions that require *Evaluation* or *Creation*.",
            llm_config=llm_config
)

case_study = autogen.AssistantAgent(
            name="case_study",
            is_termination_msg=termination_msg,
            system_message="""You are a Case Study Generation Agent. Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 15-mark case study question, create a concise and structured case study that includes:

            Back-Ground: A brief background of the scenario.
            Problem Statement: A clear description of the challenge.
            Supporting Data: Relevant information, diagrams, or examples.
            Questions: Thought-provoking questions that require learners to analyze, evaluate, and create solutions based on the scenario.(15 marks - Total)
            Ensure the case study is practical, aligned with the context, and encourages critical thinking.""",
            human_input_mode="ALWAYS",

             llm_config=llm_config,
        )

def get_problem():
        """Get the problem statement."""
        return """From the provided document, extract the most relevant topics and subtopics. 
        For each topic, retrieve the corresponding explanations, definitions, or related content. 
        Focus on identifying both *factual* and *conceptual* content that could span across 
        different cognitive levels (e.g., factual recall, application, analysis).
        Generate a question paper based on the provided document.

        If Mathematical: Extract equations and concepts to create short-answer, long-answer, 
        and problem-solving questions, ensuring balanced difficulty and a total of 100 marks.
        If Non-Mathematical: Identify key topics, definitions, and concepts to create a mix of 
        factual, application, and analysis questions, including diagrams and scenarios, 
        with a total of 100 marks.
        Ensure the format is clear and questions are categorized by type and difficulty."""

# Set up group chat
groupchat = autogen.GroupChat(
    agents=[user_proxy, ragproxyagent, short_answer, long_answer, case_study, final_paper],
    messages=[],
    max_round=9,
    speaker_transitions_type="allowed",
    speaker_selection_method="auto",
)

manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=llm_config)

manager.reset()
     
output = ragproxyagent.initiate_chat(
    manager,
    problem=get_problem(),
    message=ragproxyagent.message_generator,
    n_results=3,
)


# Define avatars for chat interface
avatar = {
    user_proxy.name: "👨‍💼",
    ragproxyagent.name: "👩‍💻",
    short_answer.name: "👩‍🔬",
    long_answer.name: "🗓",
    case_study.name: "🛠",
    final_paper.name: '📝'
}

def print_messages(recipient, messages, sender, config):
    content = messages[-1]['content']
    if all(key in messages[-1] for key in ['name']):
        chat_interface.send(content, user=messages[-1]['name'], avatar=avatar[messages[-1]['name']], respond=False)
    else:
        chat_interface.send(content, user=recipient.name, avatar=avatar[recipient.name], respond=False)
    return False, None

# Register reply functions for all agents
for agent in [user_proxy, ragproxyagent, short_answer,long_answer, case_study, final_paper]:
    agent.register_reply(
        [autogen.Agent, None],
        reply_func=print_messages,
        config={"callback": None},
    )

# Set up Panel chat interface
pn.extension(design="material")
initiate_chat_task_created = False

async def delayed_initiate_chat(agent, recipient, message):
    global initiate_chat_task_created
    initiate_chat_task_created = True
    await asyncio.sleep(2)
    await agent.a_initiate_chat(recipient, message=message)

async def callback(contents: str, user: str, instance: pn.chat.ChatInterface):
    global initiate_chat_task_created
    global input_future
    
    if not initiate_chat_task_created:
        asyncio.create_task(delayed_initiate_chat(user_proxy, manager, contents))
    else:
        if input_future and not input_future.done():
            input_future.set_result(contents)
        else:
            print("There is currently no input being awaited.")

chat_interface = pn.chat.ChatInterface(callback=callback)
chat_interface.send("Send a message!", user="System", respond=False)
chat_interface.servable()
