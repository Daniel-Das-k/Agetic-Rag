from autogen import AssistantAgent, UserProxyAgent, config_list_from_json, GroupChat, GroupChatManager
import autogen


gemini_config_list = config_list_from_json(
    "OAI_CONFIG_LIST.json",
    filter_dict={
        "model": ["gemini-2.0-flash-exp"],
    },
)

llm_config = {
    "config_list": gemini_config_list,
    "seed": 53,
    "temperature": 0.7,
    "timeout": 300
}
user_proxy = UserProxyAgent(
    name="Admin",
    system_message="A human admin. Provide a topic summary to generate Bloom's Taxonomy-based questions.",
    code_execution_config={"work_dir": "workspace", "use_docker": False},
    human_input_mode="NEVER",
    default_auto_reply="Waiting for your summary..."
)

question_generator = AssistantAgent(
    name="QuestionGenerator",
    llm_config=llm_config,
    system_message='''Question Generator. Generate questions based on Bloom's Taxonomy levels for a given topic 
    summary. Follow this structure for question generation:

    1. **Remember**: Basic recall questions (e.g., Define, List, State)
    2. **Understand**: Comprehension questions (e.g., Explain, Describe, Discuss)
    3. **Apply**: Practical application questions (e.g., Solve, Use, Demonstrate)
    4. **Analyze**: Critical thinking questions (e.g., Compare, Contrast, Differentiate)
    5. **Evaluate**: Judgement-based questions (e.g., Justify, Critique, Assess)
    6. **Create**: Synthesis questions requiring original work (e.g., Design, Formulate, Construct)

    Ensure the questions align with the provided topic and are educationally effective.'''
)


critic = AssistantAgent(
    name="Critic",
    llm_config=llm_config,
    system_message="Critic. Review the generated questions for alignment with Bloom's Taxonomy levels and accuracy. Provide constructive feedback."
)


executor = AssistantAgent(
    name="Executor",
    system_message="Executor. Compile the generated questions into a structured format and verify usability.",
    code_execution_config={"work_dir": "executor_workspace", "use_docker": False},
    llm_config=llm_config
)


planner = AssistantAgent(
    name="Planner",
    system_message='''Planner. Plan the workflow for generating Bloom's Taxonomy-based questions from a topic 
    summary. Workflow steps:

    1. Parse the topic summary.
    2. Instruct the QuestionGenerator to create questions for each Bloom's Taxonomy level.
    3. Request the Critic to validate the questions.
    4. Forward validated questions to the Executor for compilation and delivery.''',
    llm_config=llm_config
)



groupchat = GroupChat(
  
    agents=[user_proxy, question_generator, critic, executor, planner], messages=[]
)
manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)

user_proxy.initiate_chat(manager, message="The topic summary is 'Neural Networks are machine learning models inspired by the structure and function of biological neural networks. They consist of layers of nodes, which transform input data using weighted connections and activation functions. Applications include image recognition, language processing, and autonomous systems.' Generate Bloom's Taxonomy questions.")



# from autogen import AssistantAgent, UserProxyAgent, config_list_from_json, GroupChat, GroupChatManager
# import autogen

# gemini_config_list = config_list_from_json(
#     "OAI_CONFIG_LIST.json",
#     filter_dict={"model": ["gemini-2.0-flash-exp"]},
# )

# llm_config = {
#     "config_list": gemini_config_list,
#     "seed": 53,
#     "temperature": 0,
#     "timeout": 300,
# }

# user_proxy = UserProxyAgent(
#     name="Admin",
#     system_message="A human admin. Provide a topic summary to generate Bloom's Taxonomy-based questions.",
#     human_input_mode="NEVER",
#     default_auto_reply="Provide your topic summary...",
# )


# question_generator = AssistantAgent(
#     name="QuestionGenerator",
#     llm_config=llm_config,
#     system_message='''Question Generator. Generate questions based on Bloom's Taxonomy levels for a given topic 
#     summary. Follow this structure:

#     1. **Remember**: Basic recall questions (e.g., Define, List, State)
#     2. **Understand**: Comprehension questions (e.g., Explain, Describe, Discuss)
#     3. **Apply**: Practical application questions (e.g., Solve, Use, Demonstrate)
#     4. **Analyze**: Critical thinking questions (e.g., Compare, Contrast, Differentiate)
#     5. **Evaluate**: Judgement-based questions (e.g., Justify, Critique, Assess)
#     6. **Create**: Synthesis questions requiring original work (e.g., Design, Formulate, Construct)
#     '''
# )

# critic = AssistantAgent(
#     name="Critic",
#     llm_config=llm_config,
#     system_message="Critic. Validate the generated questions for accuracy and alignment with Bloom's Taxonomy.",
# )


# executor = AssistantAgent(
#     name="Executor",
#     llm_config=llm_config,
#     system_message="Executor. Compile the validated questions into a structured format for final delivery.",
# )


# planner = AssistantAgent(
#     name="Planner",
#     llm_config=llm_config,
#     system_message='''Planner. Plan the workflow to generate Bloom's Taxonomy-based questions. Steps:
#     1. Parse the topic summary.
#     2. Instruct QuestionGenerator to create questions.
#     3. Request Critic to validate questions.
#     4. Forward validated questions to Executor for final compilation.
#     ''',
# )


# groupchat = GroupChat(
#     agents=[user_proxy, question_generator, critic, executor, planner],
#     messages=[]
# )
# manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)


# user_proxy.initiate_chat(
#     manager,
#     message=""" The topic summary is " can feed the scan into a 3D printing machine and after fine-tuning the materials and the ear shapes, print the entire hearing aids." Generate Bloom's Taxonomy questions.""",
# )
