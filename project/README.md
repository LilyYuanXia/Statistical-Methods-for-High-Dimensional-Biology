# Group Project - individual report: 

Survival and Gene Expressions Analysis of Pancreatic adenocarcinoma (PAAD)

## A concise summary of contributions of each group member (1 pts)

Group tasks were divided into five steps, including project topic background information research, Data Dexcrption and Wrangling, Exploratory Data Analysis (EDA), Modeliing Building, and results interpretation. And our group members include Hassan Ali, Simran Samra, Sophia Li and me. Hassian and Simran are from Biology background, and one works in BC cancer and another has experience in Heart Lung Innovation. Sophia and I are from Statistical background. Both of us have strong data analysis coding stills but lack of the experience in analyze high dimensional Genomic data. 

Based on the each member's academic backgroup, Hassan and Simran are assigned to the background information research ("blue task"), on the other hand, Sophia and I mainly focus on the model building ("yellow task"). Data Dexcrption and Wrangling, EDA and results interpretation ("green tasks") were done all together. 

![group tasks](https://github.com/STAT540-UBC/zz_Xia-Lily_STAT540_2020/blob/master/project/plot/group%20tasks.png)

"blue task": 

1. research the backgroud of Pancreatic adenocarcinoma (PAAD) and introduce the topic to other group members
2. evaluate whether group objectives are doable or not

"green task":

1. Data Description and Wrangling: decribe the data variables (gene data and metadata) and sanitize the data
2. EDA: explore the data without any assumption. For instance, using density plots, boxplots and principal components analysis (PCA) to visualize the data and compare different cohort of subjects. 
3. Results interpratation: explain the EDA and statistical models outcome and summary analysis limitation

"yellow task":

1. survival analysis
2. vital status classification
3. Hypothesis test
4. Clustering samples and genes

## Your specific contributions and comments (1 pts)
As mention above, I mainly focus on data cleaning and data analyzing. Besides, I also help to maintain the Github repo. Overall, all tasks worked well but final results interpration and limitation discussion were not as what I expected. From my perspective, having members from genome or biology backgroud could help us not only interprate the results statistically, but also assist us to think the results in different aspects. For example, I know that although our models shows there is no statistical significant between patients' gene expression and vital status, it does not mean there is no association between them. Therefore, In addition to highly expressed genes, I hope I could have more information on what mutational genes we should mainly focus on based on backgroud research. This project is really good hand-on experience to apply the approaches that we have learnt in class and seminars, and it's fun to work with different background students. 

## Answer to one question specific to your project (3 points)
What is your rationale for mostly analyzing the top 10 genes, and not for instance genes that pass FDR, or some other number of top genes? How do you think this impacted your results?

