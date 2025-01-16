from autogen import AssistantAgent, UserProxyAgent, config_list_from_json, GroupChat, GroupChatManager
import autogen

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

user_proxy = UserProxyAgent(
    name="Admin",
    system_message="A human admin. Provide a topic summary to generate Bloom's Taxonomy-based questions.",
    human_input_mode="NEVER",
    default_auto_reply="Provide your topic summary...",
)


question_generator = AssistantAgent(
    name="QuestionGenerator",
    llm_config=llm_config,
    system_message='''Question Generator. Generate questions based on Bloom's Taxonomy levels for a given topic 
    summary. Follow this structure:

    1. **Remember**: Basic recall questions (e.g., Define, List, State)
    2. **Understand**: Comprehension questions (e.g., Explain, Describe, Discuss)
    3. **Apply**: Practical application questions (e.g., Solve, Use, Demonstrate)
    4. **Analyze**: Critical thinking questions (e.g., Compare, Contrast, Differentiate)
    5. **Evaluate**: Judgement-based questions (e.g., Justify, Critique, Assess)
    6. **Create**: Synthesis questions requiring original work (e.g., Design, Formulate, Construct)
    '''
)

critic = AssistantAgent(
    name="Critic",
    llm_config=llm_config,
    system_message="Critic. Validate the generated questions for accuracy and alignment with Bloom's Taxonomy.",
)


executor = AssistantAgent(
    name="Executor",
    llm_config=llm_config,
    system_message="Executor. Compile the validated questions into a structured format for final delivery.",
)


planner = AssistantAgent(
    name="Planner",
    llm_config=llm_config,
    system_message='''Planner. Plan the workflow to generate Bloom's Taxonomy-based questions. Steps:
    1. Parse the topic summary.
    2. Instruct QuestionGenerator to create questions.
    3. Request Critic to validate questions.
    4. Forward validated questions to Executor for final compilation.
    ''',
)


groupchat = GroupChat(
    agents=[user_proxy, question_generator, critic, executor, planner],
    messages=[]
)
manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)


user_proxy.initiate_chat(
    manager,
    message=""" The topic summary is "Linear regression models are often fitted using the least squares approach, but they may also be fitted in other ways, such as by minimizing the "lack of fit" in some other norm (as with least absolute deviations regression), or by minimizing a penalized version of the least squares cost function as in ridge regression (L2-norm penalty) and lasso (L1-norm penalty). Use of the Mean Squared Error (MSE) as the cost on a dataset that has many large outliers, can result in a model that fits the outliers more than the true data due to the higher importance assigned by MSE to large errors. So, cost functions that are robust to outliers should be used if the dataset has many large outliers. Conversely, the least squares approach can be used to fit models that are not linear models. Thus, although the terms "least squares" and "linear model" are closely linked, they are not synonymous.

Formulation

In linear regression, the observations (red) are assumed to be the result of random deviations (green) from an underlying relationship (blue) between a dependent variable (y) and an independent variable (x).
Given a data set 
{
y
i
,
x
i
1
,
…
,
x
i
p
}
i
=
1
n
{\displaystyle \{y_{i},\,x_{i1},\ldots ,x_{ip}\}_{i=1}^{n}} of n statistical units, a linear regression model assumes that the relationship between the dependent variable y and the vector of regressors x is linear. This relationship is modeled through a disturbance term or error variable ε—an unobserved random variable that adds "noise" to the linear relationship between the dependent variable and regressors. Thus the model takes the form
y
i
=
β
0
+
β
1
x
i
1
+
⋯
+
β
p
x
i
p
+
ε
i
=
x
i
T
β
+
ε
i
,
i
=
1
,
…
,
n
,
{\displaystyle y_{i}=\beta {0}+\beta _{1}x{i1}+\cdots +\beta {p}x{ip}+\varepsilon _{i}=\mathbf {x} _{i}^{\mathsf {T}}{\boldsymbol {\beta }}+\varepsilon _{i},\qquad i=1,\ldots ,n,}where T denotes the transpose, so that xiTβ is the inner product between vectors xi and β.

Often these n equations are stacked together and written in matrix notation as

y
=
X
β
+
ε
,
{\displaystyle \mathbf {y} =\mathbf {X} {\boldsymbol {\beta }}+{\boldsymbol {\varepsilon }},\,}
where

y
=
[
y
1
y
2
⋮	
y
n
 
]
,
{\displaystyle \mathbf {y} ={\begin{bmatrix}y_{1}\\y_{2}\\\vdots \\y_{n}\end{bmatrix}},\quad }
X
=
[
x
1
T
x
2
T
⋮	
x
n
T
 
]
=
[
1	
x
11
⋯	
x
1
p
1	
x
21
⋯	
x
2
p
⋮	⋮	⋱	⋮	1	
x
n
1
⋯	
x
n
p
 
]
,
{\displaystyle \mathbf {X} ={\begin{bmatrix}\mathbf {x} {1}^{\mathsf {T}}\\\mathbf {x} _{2}^{\mathsf {T}}\\\vdots \\\mathbf {x} _{n}^{\mathsf {T}}\end{bmatrix}}={\begin{bmatrix}1&x{11}&\cdots &x_{1p}\\1&x_{21}&\cdots &x_{2p}\\\vdots &\vdots &\ddots &\vdots \\1&x_{n1}&\cdots &x_{np}\end{bmatrix}},}
β
=
[
β
0
β
1
β
2
⋮	
β
p
 
]
,
ε
=
[
ε
1
ε
2
⋮	
ε
n
 
]
.
{\displaystyle {\boldsymbol {\beta }}={\begin{bmatrix}\beta _{0}\\\beta _{1}\\\beta _{2}\\\vdots \\\beta _{p}\end{bmatrix}},\quad {\boldsymbol {\varepsilon }}={\begin{bmatrix}\varepsilon _{1}\\\varepsilon _{2}\\\vdots \\\varepsilon _{n}\end{bmatrix}}.}
Notation and terminology
y
{\displaystyle \mathbf {y} } is a vector of observed values 
y
i
 
(
i
=
1
,
…
,
n
)
{\displaystyle y_{i}\ (i=1,\ldots ,n)} of the variable called the regressand, endogenous variable, response variable, target variable, measured variable, criterion variable, or dependent variable. This variable is also sometimes known as the predicted variable, but this should not be confused with predicted values, which are denoted 
y
^
{\displaystyle {\hat {y}}}. The decision as to which variable in a data set is modeled as the dependent variable and which are modeled as the independent variables may be based on a presumption that the value of one of the variables is caused by, or directly influenced by the other variables. Alternatively, there may be an operational reason to model one of the variables in terms of the others, in which case there need be no presumption of causality.
X
{\displaystyle \mathbf {X} } may be seen as a matrix of row-vectors 
x
i
⋅
{\displaystyle \mathbf {x} _{i\cdot }} or of n-dimensional column-vectors 
x
⋅
j
{\displaystyle \mathbf {x} _{\cdot j}}, which are known as regressors, exogenous variables, explanatory variables, covariates, input variables, predictor variables, or independent variables (not to be confused with the concept of independent random variables). The matrix 
X
{\displaystyle \mathbf {X} } is sometimes called the design matrix.
Usually a constant is included as one of the regressors. In particular, 
x
i
0
=
1
{\displaystyle x_{i0}=1} for 
i
=
1
,
…
,
n
{\displaystyle i=1,\ldots ,n}. The corresponding element of β is called the intercept. Many statistical inference procedures for linear models require an intercept to be present, so it is often included even if theoretical considerations suggest that its value should be zero.
Sometimes one of the regressors can be a non-linear function of another regressor or of the data values, as in polynomial regression and segmented regression. The model remains linear as long as it is linear in the parameter vector β.
The values xij may be viewed as either observed values of random variables Xj or as fixed values chosen prior to observing the dependent variable. Both interpretations may be appropriate in different cases, and they generally lead to the same estimation procedures; however different approaches to asymptotic analysis are used in these two situations.
β
{\displaystyle {\boldsymbol {\beta }}} is a 
(
p
+
1
)
{\displaystyle (p+1)}-dimensional parameter vector, where 
β
0
{\displaystyle \beta _{0}} is the intercept term (if one is included in the model—otherwise 
β
{\displaystyle {\boldsymbol {\beta }}} is p-dimensional). Its elements are known as effects or regression coefficients (although the latter term is sometimes reserved for the estimated effects). In simple linear regression, p=1, and the coefficient is known as regression slope. Statistical estimation and inference in linear regression focuses on β. The elements of this parameter vector are interpreted as the partial derivatives of the dependent variable with respect to the various independent variables.
ε
{\displaystyle {\boldsymbol {\varepsilon }}} is a vector of values 
ε
i
{\displaystyle \varepsilon _{i}}. This part of the model is called the error term, disturbance term, or sometimes noise (in contrast with the "signal" provided by the rest of the model). This variable captures all other factors which influence the dependent variable y other than the regressors x. The relationship between the error term and the regressors, for example their correlation, is a crucial consideration in formulating a linear regression model, as it will determine the appropriate estimation method.
Fitting a linear model to a given data set usually requires estimating the regression coefficients 
β
{\displaystyle {\boldsymbol {\beta }}} such that the error term
ε
=
y
−
X
β
{\displaystyle {\boldsymbol {\varepsilon }}=\mathbf {y} -\mathbf {X} {\boldsymbol {\beta }}} is minimized. For example, it is common to use the sum of squared errors 
‖
ε
‖
2
2
{\displaystyle \|{\boldsymbol {\varepsilon }}\|_{2}^{2}} as a measure of 
ε
{\displaystyle {\boldsymbol {\varepsilon }}} for minimization.

Example
Consider a situation where a small ball is being tossed up in the air and then we measure its heights of ascent hi at various moments in time ti. Physics tells us that, ignoring the drag, the relationship can be modeled as

h
i
=
β
1
t
i
+
β
2
t
i
2
+
ε
i
,
{\displaystyle h_{i}=\beta {1}t{i}+\beta {2}t{i}^{2}+\varepsilon _{i},}
where β1 determines the initial velocity of the ball, β2 is proportional to the standard gravity, and εi is due to measurement errors. Linear regression can be used to estimate the values of β1 and β2 from the measured data. This model is non-linear in the time variable, but it is linear in the parameters β1 and β2; if we take regressors xi = (xi1, xi2)  = (ti, ti2), the model takes on the standard form

h
i
=
x
i
T
β
+
ε
i
.
{\displaystyle h_{i}=\mathbf {x} _{i}^{\mathsf {T}}{\boldsymbol {\beta }}+\varepsilon _{i}.}
Assumptions
See also: Ordinary least squares § Assumptions
Standard linear regression models with standard estimation techniques make a number of assumptions about the predictor variables, the response variable and their relationship. Numerous extensions have been developed that allow each of these assumptions to be relaxed (i.e. reduced to a weaker form), and in some cases eliminated entirely. Generally these extensions make the estimation procedure more complex and time-consuming, and may also require more data in order to produce an equally precise model.[citation needed]


Example of a cubic polynomial regression, which is a type of linear regression. Although polynomial regression fits a curve model to the data, as a statistical estimation problem it is linear, in the sense that the regression function E(y | x) is linear in the unknown parameters that are estimated from the data. For this reason, polynomial regression is considered to be a special case of multiple linear regression.
The following are the major assumptions made by standard linear regression models with standard estimation techniques (e.g. ordinary least squares):

Weak exogeneity. This essentially means that the predictor variables x can be treated as fixed values, rather than random variables. This means, for example, that the predictor variables are assumed to be error-free—that is, not contaminated with measurement errors. Although this assumption is not realistic in many settings, dropping it leads to significantly more difficult errors-in-variables models.
Linearity. This means that the mean of the response variable is a linear combination of the parameters (regression coefficients) and the predictor variables. Note that this assumption is much less restrictive than it may at first seem. Because the predictor variables are treated as fixed values (see above), linearity is really only a restriction on the parameters. The predictor variables themselves can be arbitrarily transformed, and in fact multiple copies of the same underlying predictor variable can be added, each one transformed differently. This technique is used, for example, in polynomial regression, which uses linear regression to fit the response variable as an arbitrary polynomial function (up to a given degree) of a predictor variable. With this much flexibility, models such as polynomial regression often have "too much power", in that they tend to overfit the data. As a result, some kind of regularization must typically be used to prevent unreasonable solutions coming out of the estimation process. Common examples are ridge regression and lasso regression. Bayesian linear regression can also be used, which by its nature is more or less immune to the problem of overfitting. (In fact, ridge regression and lasso regression can both be viewed as special cases of Bayesian linear regression, with particular types of prior distributions placed on the regression coefficients.)

Visualization of heteroscedasticity in a scatter plot against 100 random fitted values using Matlab
Constant variance (a.k.a. homoscedasticity). This means that the variance of the errors does not depend on the values of the predictor variables. Thus the variability of the responses for given fixed values of the predictors is the same regardless of how large or small the responses are. This is often not the case, as a variable whose mean is large will typically have a greater variance than one whose mean is small. For example, a person whose income is predicted to be $100,000 may easily have an actual income of $80,000 or $120,000—i.e., a standard deviation of around $20,000—while another person with a predicted income of $10,000 is unlikely to have the same $20,000 standard deviation, since that would imply their actual income could vary anywhere between −$10,000 and $30,000. (In fact, as this shows, in many cases—often the same cases where the assumption of normally distributed errors fails—the variance or standard deviation should be predicted to be proportional to the mean, rather than constant.) The absence of homoscedasticity is called heteroscedasticity. In order to check this assumption, a plot of residuals versus predicted values (or the values of each individual predictor) can be examined for a "fanning effect" (i.e., increasing or decreasing vertical spread as one moves left to right on the plot). A plot of the absolute or squared residuals versus the predicted values (or each predictor) can also be examined for a trend or curvature. Formal tests can also be used; see Heteroscedasticity. The presence of heteroscedasticity will result in an overall "average" estimate of variance being used instead of one that takes into account the true variance structure. This leads to less precise (but in the case of ordinary least squares, not biased) parameter estimates and biased standard errors, resulting in misleading tests and interval estimates. The mean squared error for the model will also be wrong. Various estimation techniques including weighted least squares and the use of heteroscedasticity-consistent standard errors can handle heteroscedasticity in a quite general way. Bayesian linear regression techniques can also be used when the variance is assumed to be a function of the mean. It is also possible in some cases to fix the problem by applying a transformation to the response variable (e.g., fitting the logarithm of the response variable using a linear regression model, which implies that the response variable itself has a log-normal distribution rather than a normal distribution)." Generate Bloom's Taxonomy questions.""",
)
