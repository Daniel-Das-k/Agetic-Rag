# import chromadb
# import os
# import autogen
# from autogen import AssistantAgent
# from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
# from autogen import config_list_from_json
# from langchain_google_genai import GoogleGenerativeAIEmbeddings
# from langchain_chroma import Chroma
# from flask import Flask, request, jsonify,render_template
# import os
# import io
# from reportlab.lib.pagesizes import letter
# from reportlab.pdfgen import canvas
# from reportlab.lib import colors
# from textwrap import wrap

# from langchain.text_splitter import RecursiveCharacterTextSplitter

# text_splitter = RecursiveCharacterTextSplitter(separators=["\n\n", "\n", "\r", "\t"])


# app = Flask(__name__)

# class CustomEmbeddingFunction:
#     def __init__(self, embedding_model):
#         self.embedding_model = embedding_model

#     def __call__(self, input: list) -> list:
    
#         return self.embedding_model.embed_documents(input)

# gemini_config_list = config_list_from_json(
#     "OAI_CONFIG_LIST.json",
#     filter_dict={"model": ["gemini-2.0-flash-exp"]},
# )

# llm_config = {
#     "config_list": gemini_config_list,
#     "cache_seed":None,
#     "temperature": 0,
#     "timeout": 300,
# }

# def termination_msg(x):
#     return isinstance(x, dict) and "TERMINATE" == str(x.get("content", ""))[-9:].upper()

# Final_paper = AssistantAgent(
#     name="Question_format",
#     system_message="""Organize the following questions into a finalized question paper. The paper should include:
# 1. *Short-Answer Questions (2 Marks)* focusing on Remember and Understand.
# 2. *Long-Answer Questions (10 Marks)* focusing on Evaluate and Create.
# 3. *Diagram & Scenario-Based Questions (10 Marks)* focusing on Apply and Create.
# Ensure the questions are categorized, formatted appropriately, and distributed according to Bloom’s Taxonomy, with a logical structure.
# """,
#     llm_config=llm_config
# )

# short_answer = autogen.AssistantAgent(
#     name="2_marks",
#     is_termination_msg=termination_msg,
#     system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 2-mark short-answer questions focusing on *Remembering* (recalling facts) and *Understanding* (explaining concepts). These questions should test knowledge recall and comprehension, ensuring answers are brief (2-3 sentences).",
#     llm_config=llm_config,
# )

# long_answer = autogen.AssistantAgent(
#     name="10_marks",
#     is_termination_msg=termination_msg,
#     system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 10-mark long-answer questions that require *Evaluation* (critically assessing ideas or arguments) or *Creation* (designing new ideas, developing models). Ensure these questions are multi-part and require a detailed response.",
#     llm_config=llm_config,
# )

# case_study = autogen.AssistantAgent(
#     name="case_study",
#     is_termination_msg=termination_msg,
#     system_message="""You are a Case Study Generation Agent. Based on the given context, create a concise and structured case study that includes:

#             Introduction: A brief background of the scenario.
#             Problem Statement: A clear description of the challenge.
#             Supporting Data: Relevant information, diagrams, or examples.
#             Questions: Thought-provoking questions that require learners to analyze, evaluate, and create solutions based on the scenario.
#             Ensure the case study is practical, aligned with the context, and encourages critical thinking.""",
#     llm_config=llm_config,
# )


# CHROMA_DB_PATH = "/tmp/db"
# CHROMA_COLLECTION = "autogen-docs-test"

# chroma_client = chromadb.PersistentClient(path=CHROMA_DB_PATH)
# collection = chroma_client.get_or_create_collection(name=CHROMA_COLLECTION)

# embeddings = GoogleGenerativeAIEmbeddings(
#     model="models/embedding-001",
#     google_api_key="AIzaSyBt3JTEUjrNFb392jP-2YExCCSOlv8ev7A"
# )
# custom_embedding = CustomEmbeddingFunction(embedding_model=embeddings)

# vector_db = Chroma(embedding_function=custom_embedding)

# ragproxyagent = RetrieveUserProxyAgent(
#     name="ragproxyagent",
#     human_input_mode="NEVER",
#     llm_config=llm_config,
#     code_execution_config=False,
#     retrieve_config={
#         "model": gemini_config_list[0]["model"],
#         "task": "qa",
#         "update_context": True,
#         "n_results": 3,
#         "docs_path": [
#             "history.pdf"           
#         ],
#         "chunk_token_size": 3000,
#         "chunk_mode": "one_line",
#         "client": chromadb.PersistentClient(path=CHROMA_DB_PATH),
#         "get_or_create": True,
#         "overwrite": False,
#         "vector_db": vector_db,
#         "collection_name": CHROMA_COLLECTION,
#         "embedding_function": custom_embedding,  
#         "custom_text_split_function": text_splitter.split_text,
#     },
# )

# # PROBLEM = "Generate mathematical questions related to equations the question must be based on equaltions. The overall question paper should have a mix of short-answer, long-answer, and diagram-based questions. The questions should be categorized based on the level of difficulty and complexity, ensuring a balanced distribution of marks.The question paper should have max marks : 100 marks"
# PROBLEM = """From the provided document, extract the most relevant topics and subtopics. For each topic, retrieve the corresponding explanations, definitions, or related content. Focus on identifying both *factual* and *conceptual* content that could span across different cognitive levels (e.g., factual recall, application, analysis).
# Generate a question paper based on the provided document.

# If Mathematical: Extract equations and concepts to create short-answer, long-answer, and problem-solving questions, ensuring balanced difficulty and a total of 100 marks.
# If Non-Mathematical: Identify key topics, definitions, and concepts to create a mix of factual, application, and analysis questions, including diagrams and scenarios, with a total of 100 marks.
# Ensure the format is clear and questions are categorized by type and difficulty.
# """

# groupchat = autogen.GroupChat(
#     agents=[ragproxyagent, short_answer, long_answer, case_study,Final_paper], messages=[], max_round=9, speaker_selection_method="round_robin"
# )
# manager = autogen.GroupChatManager(groupchat=groupchat, llm_config=llm_config)

# # @app.route('/run', methods=['GET', 'POST'])
# # def run():
# #     print("Running")
# #     file = request.files.get('file')
# #     print(file)
# #     output = ragproxyagent.initiate_chat(
# #         manager,
# #         problem=PROBLEM,
# #         message=ragproxyagent.message_generator,
# #         n_results=3,
# #     )

# #     chat_history = output.chat_history

# #     final_output = None
# #     for message in chat_history:
# #         if 'content' in message:
# #             final_output = message['content']

# #     print(final_output)
# #     return jsonify({"final_output": final_output})


# @app.route('/run', methods=['GET', 'POST'])
# def run():
#     manager.reset()
#     print("Running")
#     file = request.files.get('file')
#     print(file)

#     output = ragproxyagent.initiate_chat(
#         manager,
#         problem=PROBLEM,
#         message=ragproxyagent.message_generator,
#         n_results=3,
#     )

#     chat_history = output.chat_history

#     final_output = None
#     for message in chat_history:
#         if 'content' in message:
#             final_output = message['content']

#     print(final_output)

#     pdf_path = os.path.join("db", "output.pdf")  

#     pdf_buffer = io.BytesIO()  

#     c = canvas.Canvas(pdf_buffer, pagesize=letter)  
#     c.setFont("Helvetica", 10) 
#     width, height = letter  

#     c.setFont("Helvetica-Bold", 12)
#     c.drawString(100, height - 40, "Generated Question Paper")
    
#     c.setFont("Helvetica", 10)  

#     wrapped_text = wrap(final_output, width=90) 
    
#     y_position = height - 60  
    
#     for line in wrapped_text:
#         c.drawString(100, y_position, line) 
#         y_position -= 12 

#         if y_position < 40:
#             c.showPage() 
#             c.setFont("Helvetica", 10) 
#             y_position = height - 40 

#     c.showPage()
#     c.save()  

#     pdf_buffer.seek(0)  

#     with open(pdf_path, "wb") as f:
#         f.write(pdf_buffer.read())

#     print(f"PDF saved to: {pdf_path}")
    
#     return jsonify({"message": "PDF saved successfully", "pdf_path": pdf_path})

# @app.route('/')
# def index():
#     return render_template('index.html')

# app.run(debug=True,port=9090)



import os
import io
from textwrap import wrap
from flask import Flask, request, jsonify, render_template
import chromadb
import autogen
import uuid
from autogen import AssistantAgent
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent
from autogen import config_list_from_json
from langchain_google_genai import GoogleGenerativeAIEmbeddings
from langchain_chroma import Chroma
from langchain.text_splitter import RecursiveCharacterTextSplitter
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib import colors

app = Flask(__name__)

class CustomEmbeddingFunction:
    def __init__(self, embedding_model):
        self.embedding_model = embedding_model

    def __call__(self, input: list) -> list:
        return self.embedding_model.embed_documents(input)

class QuestionPaperGenerator:
    def __init__(self,file):

        unique_id=uuid.uuid4() 
        self.CHROMA_DB_PATH = f"/tmp/chroma.db-{unique_id}"
        self.CHROMA_COLLECTION = "autogen-docs-test"
        self.file = file
        
        # Setup configurations
        self.setup_config()
        self.setup_database()
        self.create_agents()
        self.setup_group_chat()

    def setup_config(self):
        """Setup basic configurations for the agents."""
        self.gemini_config_list = config_list_from_json(
            "OAI_CONFIG_LIST.json",
            filter_dict={"model": [os.getenv("MODEL")]},
        )
        
        self.llm_config = {
            "config_list": self.gemini_config_list,
            "cache_seed": None,
            "temperature": 0,
            "timeout": 300,
        }

    def setup_database(self):
        self.chroma_client = chromadb.PersistentClient(path=self.CHROMA_DB_PATH)
        self.collection = self.chroma_client.get_or_create_collection(name=self.CHROMA_COLLECTION)

        self.embeddings = GoogleGenerativeAIEmbeddings(
            model=os.getenv("EMBEDDING"),
            google_api_key=os.getenv("URL"),
        )
        self.custom_embedding = CustomEmbeddingFunction(embedding_model=self.embeddings)
        self.vector_db = Chroma(embedding_function=self.custom_embedding)

    def create_agents(self):
        text_splitter = RecursiveCharacterTextSplitter(separators=["\n\n", "\n", "\r", "\t"])
 
        self.ragproxyagent = RetrieveUserProxyAgent(
            name="ragproxyagent",
            human_input_mode="NEVER",
            llm_config=self.llm_config,
            code_execution_config=False,
            retrieve_config={
                "model": self.gemini_config_list[0]["model"],
                "task": "qa",
                "update_context": True,
                "n_results": 3,
                "docs_path": [self.file],
                "chunk_token_size": 3000,
                "chunk_mode": "one_line",
                "client": self.chroma_client,
                "get_or_create": True,
                "overwrite": False,
                "vector_db": self.vector_db,
                "collection_name": self.CHROMA_COLLECTION,
                "embedding_function": self.custom_embedding,
                "custom_text_split_function": text_splitter.split_text,
            },
        )

        self.final_paper = AssistantAgent(
            name="Question_format",
            system_message="""Organize the following questions into a finalized question paper. The paper should include:
            1. *Short-Answer Questions (2 Marks)* focusing on Remember and Understand.(10 Questions)
            2. *Long-Answer Questions (13 Marks)* focusing on Evaluate and Create.(5 Questions)
            3. *Case Study (15 Marks)* focusing on Apply and Create.(1 Question)
            Ensure the questions are categorized, formatted appropriately, and distributed according to Bloom's Taxonomy, with a logical structure.""",
            llm_config=self.llm_config
        )
        
        self.short_answer = autogen.AssistantAgent(
            name="2_marks",
            is_termination_msg=self.termination_msg,
            system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 2-mark short-answer questions focusing on *Remembering* and *Understanding*.",
            llm_config=self.llm_config,
        )
        
        self.long_answer = autogen.AssistantAgent(
            name="10_marks",
            is_termination_msg=self.termination_msg,
            system_message="Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 13-mark long-answer questions that require *Evaluation* or *Creation*.",
            llm_config=self.llm_config,
        )
        
        self.case_study = autogen.AssistantAgent(
            name="case_study",
            is_termination_msg=self.termination_msg,
            system_message="""You are a Case Study Generation Agent. Based on the retrieved content for the topics: [List of Retrieved Topics and Content], generate 15-mark case study question, create a concise and structured case study that includes:

            Introduction: A brief background of the scenario.
            Problem Statement: A clear description of the challenge.
            Supporting Data: Relevant information, diagrams, or examples.
            Questions: Thought-provoking questions that require learners to analyze, evaluate, and create solutions based on the scenario.
            Ensure the case study is practical, aligned with the context, and encourages critical thinking.""",

             llm_config=self.llm_config,
        )

    def setup_group_chat(self):
        self.groupchat = autogen.GroupChat(
            agents=[self.ragproxyagent, self.short_answer, self.long_answer, 
                   self.case_study, self.final_paper],
            messages=[],
            max_round=9,
            speaker_selection_method="round_robin"
        )
        self.manager = autogen.GroupChatManager(groupchat=self.groupchat, llm_config=self.llm_config)

    @staticmethod
    def termination_msg(x):
        return isinstance(x, dict) and "TERMINATE" == str(x.get("content", ""))[-9:].upper()

    def generate_paper(self):
        self.manager.reset()
     
        output = self.ragproxyagent.initiate_chat(
            self.manager,
            problem=self.get_problem(),
            message=self.ragproxyagent.message_generator,
            n_results=3,
        )
        
        final_output = None
        for message in output.chat_history:
            if 'content' in message:
                final_output = message['content']
        
        return final_output

    @staticmethod
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

    # def generate_pdf(self, content, output_path):
    #     """Generate PDF from the question paper content."""
    #     pdf_buffer = io.BytesIO()
    #     c = canvas.Canvas(pdf_buffer, pagesize=letter)
    #     width, height = letter

    #     # Set title
    #     c.setFont("Helvetica-Bold", 12)
    #     c.drawString(100, height - 40, "Generated Question Paper")
        
        # Set content
        # c.setFont("Helvetica", 10)
        # wrapped_text = wrap(content, width=90)
        # y_position = height - 60
        
        # for line in wrapped_text:
        #     c.drawString(100, y_position, line)
        #     y_position -= 12
        #     if y_position < 40:
        #         c.showPage()
        #         c.setFont("Helvetica", 10)
        #         y_position = height - 40

        # c.showPage()
        # c.save()
        # pdf_buffer.seek(0)

        # # Save PDF
        # os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # with open(output_path, "wb") as f:
        #     f.write(pdf_buffer.getvalue())


@app.route('/run', methods=['GET', 'POST'])
def run():
    try:
        file = request.files.get('file')
        if not file:
            return jsonify({"error": "No file uploaded"}), 400
            
        print("Running with file:", file.filename)

        generator = QuestionPaperGenerator(file=file.filename)

        final_output = generator.generate_paper()
        
        if not final_output:
            return jsonify({"error": "No output generated"}), 400
            
        pdf_path = os.path.join("db", "output.pdf")
        # generator.generate_pdf(final_output, pdf_path)
        
        # print(f"PDF saved to: {pdf_path}")
        
        return jsonify({
            "message": "PDF saved successfully",
            "pdf_path": pdf_path
        })
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=9090)