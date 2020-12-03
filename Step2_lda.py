from gensim import corpora, models
import gensim
import csv
import pandas as pd
import numpy as np

##################### STEP 1 DATA PREPROCESSING ########################

### read csv file, with each line represents each sample
f = open('1_rank2word.csv')
csv_f = csv.reader(f, delimiter=",")

# Create doc_set, total set for document
doc_set = []
for row in csv_f:
	doc_set.append(row)

# combines every word in every document to a single list
texts = []
for i in doc_set:
	texts.append(i)
  

### create dictionary
# turn our documents into a id <-> term dictionary
dictionary = corpora.Dictionary(texts)
dictionary.compactify()


### create corpus
# convert documents into a document-term matrix
class MyCorpus(object):
	def __iter__(self):
		for doc in doc_set:
			yield dictionary.doc2bow(doc)
			
corpus_memory_friendly = MyCorpus()
corpora.MmCorpus.serialize('corpus.mm', corpus_memory_friendly)
corpus = corpora.MmCorpus('corpus.mm')




##################### STEP 2 TRAINING LDA ########################

for knum in [3]:    #number of topics, append more number for iteratively train different model with different knum
	print(corpus)
	print(dictionary)
	num_unique_token = len(dictionary.token2id)
	num_documents=len(corpus)
	tfidf = models.TfidfModel(corpus)
	corpus_tfidf = tfidf[corpus]

	outtermprob = '2_k'+str(knum)+'termProb.txt'     #filename for term probability
	outdocprob  = '2_k'+str(knum)+'topicProb.txt'    #filename for document probability
    
	### Training LDA model
	print("Training lda with number of topic = %d ..." % knum)
	ldamodel = gensim.models.ldamodel.LdaModel(corpus = corpus, num_topics=knum, id2word = dictionary, passes=10)
    
    
	### Output term (word probability)
	# create empty dataframe
	data = np.array([np.zeros(num_unique_token)]*(knum*2)).T
	index = range(1, num_unique_token+1)
	columns = range(1, (knum*2)+1)
	termtopic = pd.DataFrame(data, index=index, columns=columns)
    
	# fill in the termtopic and save to file
	print("loading word probability into file " + outtermprob + "...")
	for i in range(knum):
		print("topic " + i + "...")		
		t = 1
		for item in ldamodel.get_topic_terms(i, topn=num_unique_token):	#return list of (word-id, probability) tuples
			termtopic.loc[ t, (i*2)+1] = dictionary[item[0]]        #retrieve word from dictionary based on word-id
			termtopic.loc[ t, (i*2)+2] = item[1]                    #write word possibility at another file
			t += 1
	termtopic.to_csv(outtermprob, index=False)				#output: probability of term in each topic
            
    
	### Output document (topic probability)
	# create empty dataframe
	data = np.array([np.zeros(num_documents)]*knum).T
	index = range(1, num_documents+1)
	columns = range(1, knum+1)
	doctopic = pd.DataFrame(data, index=index, columns=columns)

	# fill in the doctopic and save to file
	print("loading topic probability into file " + outdocprob + "...")
	for i in range(num_documents):
		for ind, prob in ldamodel.get_document_topics(corpus[i]):
			doctopic[ind+1][i+1] = prob
	doctopic.to_csv(outdocprob, index=False)				#output: probability of sample in each topic


